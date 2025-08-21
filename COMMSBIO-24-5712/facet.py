# Author: Ben Skubi
# Date: Jan. 30, 2025
# Facet is now maintained at https://pypi.org/project/amethyst-facet/

import numpy as np
from numpy.typing import NDArray
from typing import *
import time
import click
from dataclasses import dataclass
from functools import singledispatch
import polars as pl
import itertools
import glob
import parse
from concurrent.futures import ProcessPoolExecutor
import h5py
import sys
import math
from pathlib import Path
import duckdb
from warnings import warn

np.random.seed(0)

################################################
# Utility methods
#

def get_breakpoints(array: NDArray):
    """Get values of indexes N where array[N] != array[N+1]"""
    return np.append(0, np.where(array[:-1] != array[1:])[0] + 1)

def iter_breakpoints(
        array: NDArray[Any]
) -> Generator[Tuple[int, int], None, None]:
    """Iterate through tuples (b1, b2) where b1 and b2 are consecutive breakpoints in array

    For N breakpoints B1, B2, B3, ... BN, B1 = 0 and BN = None. N-1 tuples will be yielded. If tuple M is (BX, BX+1) then tuple M+1 is (BX+1, BX+2). 
    """
    array_breaks = get_breakpoints(array).tolist() + [None]
    
    for b1, b2 in zip(array_breaks[:-1], array_breaks[1:]):
        yield b1, b2

#############################################
# Aggregation methods

def aggsorted_uniform_disjoint_covering(
        opos: NDArray[np.int_], 
        ovalues: Dict[Any, NDArray[Any]],
        wsize: int,
        woffset: int,
        agg: np.ufunc,
) -> Tuple[NDArray[np.int_], NDArray[np.int_], Dict[Any, NDArray[Any]]]:
    """Aggregate values on windows that are the same size, gapless and non-overlapping.

    Parameters
    ----------
    opos : NDArray[np.int_]
        The observed positions array.
    ovalues : Dict[Any, NDArray[Any]]
        Dictionary mapping keys to arrays of values to be aggregated.
    wsize : int
        The fixed window size.
    woffset : int
        The offset to align windows.
    agg : np.ufunc
        The aggregation function (e.g., np.add).

    Returns
    -------
    wstarts_obs : NDArray[np.int_]
        The start positions of the aggregated windows with at least one observed value.
    wends_obs : NDArray[np.int_]
        The end positions of the aggregated windows with at least one observed value.
    ovalues_agg : Dict[Any, NDArray[Any]]
        The aggregated values for each key in `ovalues` for the windows with at least one observed value.

    """
    # Floor observed pos values to nearest bin offset (i.e. for offset = 0, binsize = 500, 0 to 499 are mapped to 0, 500-999 to 500, etc)
    owindows = ((opos - woffset) // wsize) * wsize + woffset

    # Get indices within binned positions where consecutive values differ, which will form the breakpoints for the reduceat operation
    reduceat_idx = get_breakpoints(owindows)

    # Aggregate observed values within the breakpoints, forming window sums
    ovalues_agg = {k: agg.reduceat(ovalues[k], reduceat_idx) for k in ovalues}

    # Compute a nonzero ("[key]_nz") column for each value
    orig_keys = list(ovalues.keys())
    for k in orig_keys:
        ovalues_agg[k + "_nz"] = (
            agg.reduceat(
                (ovalues[k] > 0).astype(int),
                reduceat_idx
            )
        )

    # Get the start and end positions for each window
    wstarts_obs = owindows[reduceat_idx]
    wends_obs = wstarts_obs + wsize
    return wstarts_obs, wends_obs, ovalues_agg


def aggsorted_chromwindows(
        opos: NDArray[np.int_],
        ovalues: Dict[Any, NDArray[Any]],
        wstarts: NDArray[np.int_],
        wends: NDArray[np.int_],
        agg: np.ufunc
) -> Tuple[NDArray[np.int_], NDArray[np.int_], Dict[Any, NDArray[Any]]]:
    """Compute aggregations for nonuniform chrom-specific windows (i.e. using ChromWindowsPlan)
    """
    # Store aggregated observed values
    ovalues_agg = {k: [] for k in ovalues}

    # Get the indices in the observed positions where the window start and end positions map to
    wstart_idx = np.searchsorted(opos, wstarts)
    wend_idx = np.searchsorted(opos, wends)

    # Compute window aggregations for each value
    orig_keys = ovalues_agg.keys()
    for k in orig_keys:
        ovalues_binary = (ovalues[k] > 0).astype(ovalues[k].dtype)
        # s and e store the first and 1 beyond the last index respectively of the bp-resolution values that should be aggregated for the current window
        for s, e in zip(wstart_idx, wend_idx):
            # Aggegate the values for the current window and store
            ovalues_agg[k].append(agg.reduce(ovalues[k][s:e]))

            # Aggregate the binarized values so that we sum the number of nonzero rows for the value within each window
            ovalues_agg[k + "_nz"].append(agg.reduce(ovalues_binary[k][s:e]))

    if not len(wstarts) or not len(wends) or not any(len(v) for v in ovalues_agg.values()):
        return np.array([]), np.array([]), {k: np.array([]) for k in ovalues_agg}

    nonzero = np.unique(np.concatenate([np.where(ovalues_agg[k]) for k in ovalues_agg]))
    nonzero.sort()

    return wstarts[nonzero], wends[nonzero], ovalues_agg


@dataclass
class AggregationPlan:
    agg: np.ufunc = np.add
    chrom: Optional[Any] = None
    contexts: Optional[List[str]] = None
    skipbc: Optional[List[str]] = None
    requirebc: Optional[List[str]] = None

@dataclass
class UniformDisjointCoveringPlan(AggregationPlan):
    size: int = None
    offset: int = 0


@dataclass
class UniformSteppedCoveringPlan(AggregationPlan):
    step: int = None
    size: int = None
    offset: int = 0

    def __iter__(self):
        for offset in range(self.offset, self.offset + self.size, self.step):
            yield UniformDisjointCoveringPlan(size=self.size, offset=offset, agg = self.agg)
    
    def validate(self):
        assert isinstance(self.step, int) and self.step > 0, f"{self} step size 'step' must be a positive integer"
        assert isinstance(self.size, int) and self.size > 0, f"{self} window size 'size' must be a positive integer"
        assert isinstance(self.offset, int), f"{self} offset 'offset' must be an integer"
        assert self.size % self.step == 0, f"{self} window size 'size' must be divisible by step size 'step'"

@dataclass
class ChromWindowsPlan(AggregationPlan):
    windows: Dict[Any, Tuple[NDArray, NDArray]] = None

    @property
    def chrom_wstart(self):
        if not self.chrom:
            return np.concatenate(it[0] for it in self.windows.values())
        else:
            chrom = self.chrom.decode()
            windows = self.windows.get(chrom)
            if windows:
                return windows[0]
            else:
                return np.array([])
    
    @property
    def chrom_wend(self):
        if not self.chrom:
            return np.concatenate(it[1] for it in self.windows.values())
        else:
            chrom = self.chrom.decode()
            windows = self.windows.get(chrom)
            if windows:
                return windows[1]
            else:
                return np.array([])

@singledispatch
def aggsorted(
        agg,
        opos: NDArray[np.int_],
        ovalues: Dict[Any, NDArray[Any]]
) -> Tuple[NDArray[np.int_], NDArray[np.int_], Dict[Any, NDArray[Any]]]:
    raise NotImplementedError(f"aggsorted is not implemented for agg of type {type(agg)}, but can be defined using @aggsorted.register (see functools.singledispatch)")

@aggsorted.register
def aggsorted_uniform_disjoint_covering_plan(
        plan: UniformDisjointCoveringPlan,
        opos: NDArray[np.int_],
        ovalues: Dict[Any, NDArray[Any]]
) -> Tuple[NDArray[np.int_], NDArray[np.int_], Dict[Any, NDArray[Any]]]:
    """Compute aggregations on sorted observations using uniform-size windows
    """
    return aggsorted_uniform_disjoint_covering(opos, ovalues, plan.size, plan.offset, plan.agg)

@aggsorted.register
def aggsorted_uniform_stepped_covering_plan(
        plan: UniformSteppedCoveringPlan,
        opos: NDArray[np.int_],
        ovalues: Dict[Any, NDArray[Any]] 
) -> Tuple[NDArray[np.int_], NDArray[np.int_], Dict[Any, NDArray[Any]]]:
    """Compute aggregations on sorted observations using uniform-size windows that overlap by constant step size
    """
    plan.validate()
    wstart = []
    wend = []
    wvalues = {k: [] for k in ovalues}

    for subplan in plan:
        wstart_sub, wend_sub, wvalues_sub = aggsorted(subplan, opos, ovalues)
        wstart.append(wstart_sub)
        wend.append(wend_sub)
        for k in wvalues:
            wvalues[k].append(wvalues_sub[k])
    
    wstart = np.concatenate(wstart)
    # kind='stable' seems to be about 30% faster than default kind='quicksort'
    # numpy 2.0 offers an option 'stable' (bool) which preserves relative order when values are equal. This is not necessary under the assumptions of this method because it assumes that all windows generated across all steps have distinct start sites.
    sorted_idx = np.argsort(wstart, kind='stable')
    wstart = wstart[sorted_idx]
    wend = np.concatenate(wend)[sorted_idx]
    wvalues = {k: np.concatenate(v)[sorted_idx] for k, v in wvalues.items()}

    return wstart, wend, wvalues

@aggsorted.register
def aggsorted_chromwindows_plan(
        plan: ChromWindowsPlan,
        opos: NDArray[np.int_],
        ovalues: Dict[Any, NDArray[Any]]
) -> Tuple[NDArray[np.int_], NDArray[np.int_], Dict[Any, NDArray[Any]]]:
    """Compute aggregations on sorted observations using uniform-size windows
    """
    
    return aggsorted_chromwindows(opos, ovalues, plan.chrom_wstart, plan.chrom_wend, plan.agg)    

def aggsorted_chroms(
        aggregation_plan: AggregationPlan,
        ochrom: NDArray,
        opos: NDArray[np.int_],
        ovalues: Dict[Any, NDArray[Any]],
) -> Tuple[NDArray, NDArray[np.int_], NDArray[np.int_], Dict[Any, NDArray[Any]]]:
    """Compute aggregations over a list of observations sorted by (chrom, pos)

    Returns:
        wchrom
        wstart
        wend
        wvalues 
    """
    wchrom = []
    wstart = []
    wend = []
    wvalues = {}

    # b1 and b2 iterate through indices where chromosome names change to obtain blocks of observations for a specific chromosome
    for b1, b2 in iter_breakpoints(ochrom):
        # Get the chromosome name and chrom-specific positions, and values
        chrom = ochrom[b1]
        chrom_opos = opos[b1:b2]
        chrom_ovalues = {k:ovalues[k][b1:b2] for k in ovalues}

        # Note the chromosome in the aggregation plan as some types of aggregation plans may behave differently for different chromosomes (i.e. nonuniform windows where boundary positions are chromosome-specific).
        aggregation_plan.chrom = chrom

        # Get the window start, end and values for the chromosome
        chrom_wstart, chrom_wend, chrom_wvalues = aggsorted(aggregation_plan, chrom_opos, chrom_ovalues)

        # Combine with previously observed window aggregations for other chromosomes
        if not chrom_wstart.shape[0] == 0:
            # Build a chromosome name column
            wchrom.append(np.repeat(chrom, len(chrom_wstart)))

            # Start and end positions
            wstart.append(chrom_wstart)
            wend.append(chrom_wend)

            # Put in the values. Nonzero columns suffixed with _nz are also computed in addition to the submitted values, so there are more value keys in the returned dict than the submitted values dict
            for k in chrom_wvalues:
                wvalues.setdefault(k, [])
                wvalues[k].append(chrom_wvalues[k])

    # Combine window aggregations for all chromosomes
    wchrom = np.concatenate(wchrom)
    wstart = np.concatenate(wstart)
    wend = np.concatenate(wend)
    wvalues = {k: np.concatenate(v) for k, v in wvalues.items()}

    return wchrom, wstart, wend, wvalues

###########################################################################3
#  Tests of windowing accuracy and runtime
#

def test_uniform_disjoint_covering(plan: UniformDisjointCoveringPlan, ochrom, opos, ovalues, verbose = False):
    assert plan.agg == np.add, "Tests add/sum aggregation only"

    def _aggregate(df):
        return (
            df
            .with_columns((((pl.col("pos") - plan.offset) // plan.size) * plan.size + plan.offset).alias("start"))
            .group_by("chrom", "start")
            .agg([pl.sum(k) for k in ovalues])
            .sort("chrom", "start")
            .with_columns((pl.col("start") + plan.size).alias("end"))
            .select("chrom", "start", "end", *list(ovalues.keys()))
        )

    start = time.time()
    target = _aggregate(pl.DataFrame({"chrom": ochrom, "pos": opos, **ovalues}))
    polars_time = time.time() - start

    start = time.time()
    wchrom, wstart, wend, wvalues = aggsorted_chroms(
        plan, ochrom, opos, ovalues
    )
    aggsorted_time = time.time() - start

    test = pl.DataFrame({"chrom": wchrom, "start": wstart, "end": wend, **wvalues})
    if verbose:
        print("Polars runtime (s):", polars_time, "aggsorted runtime (s):", aggsorted_time, "polars/aggsorted runtime:", polars_time/aggsorted_time, file=sys.stderr)
        print(target, file=sys.stderr)
        print(test, file=sys.stderr)
    if not test.equals(target):
        i = 0
        for target_row, test_row in zip(target.rows(), test.rows()):
            if target_row != test_row:
                print("Mismatch: Target row:", target_row, "Test row:", test_row, file=sys.stderr)
            i += 1
            if i == 10:
                break
        assert test.equals(target), f"Target and test dataframes did not match.\nTarget:\n{target}\nTest:\n{test}"

def test_uniform_stepped_covering(plan: UniformSteppedCoveringPlan, ochrom, opos, ovalues, verbose = False):
    assert plan.agg == np.add, "Tests add/sum aggregation only"

    def _aggregate(df) -> pl.DataFrame:
        dfs = []
        for subplan in plan:
            dfs.append(
                df
                .with_columns((((pl.col("pos") - subplan.offset) // subplan.size) * subplan.size + subplan.offset).alias("start"))
                .group_by("chrom", "start")
                .agg([pl.sum(k) for k in ovalues])
                .sort("chrom", "start")
                .with_columns((pl.col("start") + subplan.size).alias("end"))
                .select("chrom", "start", "end", *list(ovalues.keys()))
            )
        return (
            pl.concat(dfs)
            .sort('chrom', 'start', 'end')
        )
        

    start = time.time()
    target = _aggregate(pl.DataFrame({"chrom": ochrom, "pos": opos, **ovalues}))
    polars_time = time.time() - start

    start = time.time()
    wchrom, wstart, wend, wvalues = aggsorted_chroms(
        plan, ochrom, opos, ovalues
    )
    aggsorted_time = time.time() - start

    test = pl.DataFrame({"chrom": wchrom, "start": wstart, "end": wend, **wvalues})
    if verbose:
        print("Polars runtime (s):", polars_time, "aggsorted runtime (s):", aggsorted_time, "polars/aggsorted runtime:", polars_time/aggsorted_time, file=sys.stderr)
        print("Target:", file=sys.stderr)
        print(target, file=sys.stderr)
        print("Test:")
        print(test, file=sys.stderr)

    if not test.equals(target):
        mismatches = 0
        for i, rows in enumerate(zip(target.rows(), test.rows())):
            target_row, test_row = rows
            if target_row != test_row:
                print(f"Row {i} mismatch: Target row:", target_row, "Test row:", test_row, file=sys.stderr)
                mismatches += 1
            if mismatches == 10:
                break
        assert test.equals(target), f"Target and test dataframes did not match.\nTarget:\n{target}\nTest:\n{test}"

def test():
    ochrom = []
    opos = []
    ovalues = {k: [] for k in ['c', 't']}

    for chrom in ['chr1', 'chr2', 'chr3']:
        opos_new = np.random.randint(0, 250_000_000, 1_000_000)
        opos_new.sort()
        opos_new_unique = opos_new[get_breakpoints(opos_new)]
        opos.append(opos_new_unique)

        for k in ovalues:
            ovalues[k].append(np.random.randint(1, 3, size=len(opos_new_unique)))
        ochrom.append(np.repeat(chrom, len(opos_new_unique)))

    ochrom = np.concatenate(ochrom)
    opos = np.concatenate(opos)
    ovalues = {k: np.concatenate(v) for k, v in ovalues.items()}

    test_uniform_disjoint_covering(UniformDisjointCoveringPlan(500), ochrom, opos, ovalues)
    test_uniform_stepped_covering(UniformSteppedCoveringPlan(250, 1000), ochrom, opos, ovalues)

###############################################################
# Amethyst 
#

def aggregate_amethyst_h5_worker(
        args
):
    h5_filename, plan, origin_dataset, aggregated_dataset = args
    try:
        with h5py.File(h5_filename) as f:
            origin = f[origin_dataset]
            try:
                aggregated = aggsorted_chroms(plan, origin["chr"], origin["pos"], {'c': origin["c"], 't': origin["t"]})

                aggregated = {
                    'chr': aggregated[0],
                    'start': aggregated[1],
                    'end': aggregated[2],
                    **aggregated[3]
                }
                
            except Exception as e:
                warn(e)
                warn(f"Could not apply {plan} to {h5_filename}::{origin_dataset}:\n{origin[:10]}, writing empty dataset to output")
                
                aggregated = {
                    'chr': np.array([], dtype='S10'),
                    'start': np.array([], dtype=np.int64),
                    'end': np.array([], dtype=np.int64),
                    'c': np.array([], dtype=np.int64),
                    't': np.array([], dtype=np.int64),
                    'c_nz': np.array([], dtype=np.int64),
                    't_nz': np.array([], dtype=np.int64)
                }
    
            try:
                dtypes = [col.dtype for col in aggregated.values()]
                aggregated = np.rec.fromarrays(aggregated.values(), names=list(aggregated.keys()), formats=dtypes)
            except Exception as e:
                warn(e)
                warn(f"Could not apply {plan} to {h5_filename}::{origin_dataset}:\n{origin[:10]} (could not convert to recarray), writing empty dataset to output")
            return aggregated_dataset, aggregated
    except:
        raise ValueError(f"Could not open {h5_filename} as H5 file.")

def aggregate_amethyst_h5(
        h5_filename: str, 
        plan_name: str,
        plan: AggregationPlan,
        observations = "1",
        ignore_context = ["metadata"],
        batch_size = 1,
        overwrite: bool = True,
        compression: Optional[str] = None,
        compression_opts: Optional[Union[str,int]] = None
):
    with h5py.File(h5_filename) as f:
        datasets = {}
        for context in f:
            if context in ignore_context or plan.contexts and context not in plan.contexts:
                continue
            for barcode in f[context]:
                if (
                    (not plan.skipbc or barcode not in plan.skipbc)
                     and (not plan.requirebc or barcode in plan.requirebc)
                     and (observations in f[context][barcode])
                ):
                    barcode = f"/{context}/{barcode}/"
                    datasets[barcode + observations] = barcode + plan_name
    

    start = time.time()
    agg_time = 0
    write_time = 0
    for i, batch in enumerate(itertools.batched(datasets.items(), batch_size)):

        print(f"{h5_filename} batch {i+1}/{math.ceil(len(datasets)/batch_size)}:\n", batch, file=sys.stderr)
        batch = [(h5_filename, plan, *b) for b in batch]
        
        agg_start = time.time()
        with ProcessPoolExecutor() as ppe:
            results = list(ppe.map(aggregate_amethyst_h5_worker, batch))
        agg_time += time.time() - agg_start
        write_start = time.time()
        with h5py.File(h5_filename, 'a') as f:
                for dataset_name, data in results:
                    if dataset_name in f and overwrite:
                        del f[dataset_name]

                    try:
                        compression_opts = int(compression_opts)
                    except:
                        pass
                    f.create_dataset(dataset_name, data=data, compression = compression, compression_opts = compression_opts)
        write_time += time.time() - write_start
    runtime = time.time() - start
    print(f"Total runtime for {h5_filename} (s):", runtime, "Total aggregation time (s):", agg_time, "Total write time (s):", write_time, "Average runtime per barcode (s):", runtime/len(datasets) if len(datasets) > 0 else None, file=sys.stderr)

######################################################
# Command-line interface and CLI utility functions
#

def multiparse(formats, string, require=True):
    for format in formats:
        attempt = parse.parse(format, string)
        if attempt:
            return attempt
    if require:
        raise ValueError(f"Expected string '{string}' to parse to one of the following formats: {formats}")

def delete_from_h5(args: Tuple[str]):
    """Delete all datasets matching dataset_name for all contexts and barcodes
    """
    amethyst_h5_file, name, level = args

    def safe_del(obj):
        try:
            del obj
        except:
            pass

    with h5py.File(amethyst_h5_file, 'a') as f:
        for context in f:
            context_grp = f[f"/{context}"]

            if level == "context" and context == name:
                safe_del(context_grp)
            elif context == "metadata":
                continue
            if level in ["barcode", "dataset"]:
                for barcode in context_grp:
                    barcode_grp = context_grp[barcode]
                    if level == "barcode" and barcode == name:
                        safe_del(barcode_grp)
                    elif level == "dataset":
                        for dataset in barcode_grp:
                            if dataset == name:
                                safe_del(barcode_grp[name])

@click.group()
def facet():
    pass

@facet.command()
@click.option(
    "--compression",
    type=str,
    default = 'gzip'
)
@click.option(
    "--compression_opts",
    type=str,
    default = '6'
)
@click.argument(
    "h5_in"
)
@click.argument(
    "h5_out"
)
def convert(compression, compression_opts, h5_in, h5_out):
    """Convert an old Amethyst HDF5 file format to v2.0.0 format

    The old format stored observations as (chr, pos, pct, c, t) in a dataset at /context/barcode. The new format stores base-pair observations as (chr, pos, c, t) in a dataset at at /context/barcode/1 and window aggregations in a dataset at /context/barcode/[window_dataset_name]. It also contains a dataset at /metadata/version equal to 'amethyst2.0.0'.
    """
    try:
        compression_opts = int(compression_opts)
    except:
        pass

    if compression == "":
        compression = None
    if compression_opts == "":
        compression_opts = None

    with h5py.File(h5_in) as in_f:
        with h5py.File(h5_out, 'a') as out_f:
            # We only convert from the first format to v2.0.0
            if "/metadata/version" in in_f:
                warn("Input file contains /metadata/version, suggesting it is already converted.")
                return

            for context_str in in_f:
                context = in_f[context_str]
                for i, barcode in enumerate(context):
                    barcode_count = len(context)
                    print(f"\rConverting /{context_str}/{barcode} ({i+1}/{barcode_count})")

                    data = context[barcode][:]
                    out_f.create_dataset(f"/{context_str}/{barcode}/1", data=data, compression=compression, compression_opts=compression_opts)

            out_f.create_dataset('/metadata/version', data = "amethyst2.0.0")

@facet.command()
@click.option(
    "--h5", "--globh5", "--glob", "-g", "_globs",
    multiple=True,
    type=str,
    help = "Amethyst v2 files structured as /[context]/[barcode]/[observations]"
)
@click.option(
    "--nproc", "-p", "nproc",
    type=int,
    default=1,
    show_default = True,
    help = "Number of processes to use when aggregating multiple barcodes in a single H5 file"
)
@click.argument(
    "h5obj",
    type = click.Choice(["context", "barcode", "dataset"], case_sensitive=True)
)
@click.argument(
    "h5obj_name",
    type=str
)
@click.argument(
    "filenames",
    nargs=-1,
    type=str
)
def delete(_globs: Tuple[str], nproc: int, h5obj: str, h5obj_name: str, filenames: Tuple[str]):
    """Delete contexts, barcodes, or datasets from an Amethyst 2.0.0 format HDF5 file

Required arguments:

{context|barcode|dataset} (case sensitive), the type of object to be deleted. Skips nonexistent groups and datasets.

H5OBJ_NAME: the name of the object type to be deleted

FILENAMES: a glob or list of Amethyst v 2.0.0 filenames



Example to delete all datasets named "500" in every context and barcode:

python facet.py delete dataset 500 demo.h5
    """
    filenames: List[str] = list(filenames) + list(itertools.chain.from_iterable([glob.glob(it) for it in _globs]))
  
    with ProcessPoolExecutor(max_workers=nproc) as ppe:
        [
            ppe.submit(delete_from_h5, (filename, h5obj_name, h5obj)).result()
            for filename in filenames
        ]

@facet.command()
@click.option(
    "--h5", "--globh5", "--glob", "-g", "_globs",
    multiple=True,
    type=str,
    help = "Amethyst v2 files structured as /[context]/[barcode]/[observations]"
)
@click.option(
    "--context", "-c", "contexts",
    multiple=True,
    type=str,
    help = "Limit aggregations to these contexts. Multiple can be specified. If none given, aggregates all barcodes."
)
@click.option(
    "--skipbc",
    type = str,
    help = "Skip barcodes listed in the given newline-separated file."
)
@click.option(
    "--requirebc",
    type = str,
    help = "Require barcodes listed in the given newline-separated file."
)
@click.option(
    "--observations", "--obs", "-o", "observations",
    type=str,
    default = "1",
    show_default=True,
    help = "Name of observations dataset to aggregate in Amethyst H5 files at /[context]/[barcode]/[observations] with columns (chr, pos, c, t)"
)
@click.option(
    "--windows", "--win", "-w", "windows",
    multiple=True,
    type=str,
)
@click.option(
    "--uniform", "--unif", "--uw", "-u", "unif",
    multiple=True,
    type=str,
    help = r"Uniform window sums. Format options: {name}={size}, {name}={size}:{step}, {size}:{step}, {size} where {name} is the datasetname for the aggregation stored under /[context]/[barcode]/[name], {size} is the window size, {step} is the constant stride between window start sites (defaults to size). Window name defaults to size (if no step) or [size]by[step] (if step is specified). Examples: -u 1k=1000 -u 1k_500=1000:500 -u 1000:500 -u 1000"
)
@click.option(
    "--compression",
    type=str,
    default = 'gzip'
)
@click.option(
    "--compression_opts",
    type=str,
    default = '6'
)
@click.option(
    "--nproc", "-p", "nproc",
    type=int,
    default=1,
    show_default = True,
    help = "Number of processes to use when aggregating multiple barcodes in a single H5 file"
)
@click.argument(
    "filenames",
    nargs=-1,
    type=str
)
def agg(_globs, contexts, skipbc, requirebc, observations, unif, windows, compression, compression_opts, nproc, filenames):
    """Compute window sums over methylation observations stored in Amethyst v2 format.

    FILENAMES: Amethyst H5 filenames in format /[context]/[barcode]/[observations] to compute window sums. Can be specified as a single glob (i.e. *.h5) and will be combined with additional globs specified with -g.
    """
    filenames = list(filenames) + list(itertools.chain.from_iterable([glob.glob(it) for it in _globs]))

    plans = {}

    if skipbc:
        try:
            skipbc = [r.strip() for r in open(skipbc).readlines()]
        except:
            raise ValueError(f"Could not open or parse list of barcodes to skip at {skipbc}")
    
    if requirebc:
        try:
            requirebc = [r.strip() for r in open(requirebc).readlines()]
        except:
            raise ValueError(f"Could not open or parse list of required barcodes at {requirebc}")

    for u in unif:
        parsed = multiparse(["{name}={size}:{step}", "{name}={size}", "{size}:{step}", "{size}"], u).named
        size = parsed.get("size")
        try:
            size = int(size)
        except:
            raise ValueError(f"In {u}, could not convert size '{size}' to int")
        step = parsed.get("step", None)
        try:
            step = int(step) if step else None
        except:
            raise ValueError(f"In {u}, could not convert step '{step}' to int")
        default_name = str(size) + "_by_" + str(step) if step else str(size)
        name = parsed.get("name", default_name)
        if step:
            plan = UniformSteppedCoveringPlan(step = step, size = size, contexts=contexts, skipbc=skipbc, requirebc=requirebc)
        else:
            plan = UniformDisjointCoveringPlan(size = size, contexts=contexts, skipbc=skipbc, requirebc=requirebc)
        assert name not in plans, f"Duplicate plan {name} found: {plans[name]}, {plan}"
        plans[name] = plan

    for w in windows:
        parsed = multiparse(["{name}={filename}", "{filename}"], w).named
        filename = parsed.get("filename")
        assert filename, f"In -w {w}, filename '{filename}' could not be parsed."

        path = Path(filename)
        assert path.exists(), f"In -w {w}, filename '{filename}' does not exist."
        
        try:
            windows = duckdb.read_csv(filename).pl()
        except:
            raise ValueError(f"Unable to load windows from {filename} using DuckDB CSV sniffer.")
        
        assert set(windows.columns).intersection({"chr", "start", "end"}) == {"chr", "start", "end"}, f"For -w {w} with filename {filename}, 'chr', 'start', and 'end' columns are required but were not found -- columns were {windows.columns}.\n{windows}"

        default_filename = path.name.removesuffix(path.suffix)
        name = parsed.get("name", default_filename)

        windows = (
            windows
            .group_by("chr")
            .agg(pl.col("start"), pl.col("end"))
            .sort("chr")
        )

        windows = {
            r[0]: (np.array(r[1]), np.array(r[2]))
            for r in windows.rows()
        }

        plan = ChromWindowsPlan(windows=windows, contexts=contexts, skipbc=skipbc, requirebc=requirebc)
        plans[name] = plan

    if not plans:
        print("No plans specified.", file=sys.stderr)
    else:
        print("Plans:", file=sys.stderr)
        print(plans, file=sys.stderr)

    for filename in filenames:
        print("Aggregating", filename)
        for name, plan in plans.items():
            print(f"Writing {name}=", plan, "to", filename, file=sys.stderr)

            aggregate_amethyst_h5(
                filename, 
                name, 
                plan, 
                observations = observations, 
                batch_size=nproc, 
                compression = compression, 
                compression_opts = compression_opts
            )

            print(f"Finished writing {name}=", plan, "to", filename, file=sys.stderr)

@facet.command()
def version():
    print("Facet v. 0.1.1 (Jan. 30, 2025)")

if __name__ == "__main__":
    facet()
