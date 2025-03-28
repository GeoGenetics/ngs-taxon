import os
import pandas as pd
from typing import List, Dict
from collections import defaultdict
from snakemake.io import Namedlist


#################
### FUNCTIONS ###
#################

### General


def flatten(list_of_lists: List) -> List:
    """Flatten an irregular list of lists recursively

    https://stackoverflow.com/a/53778278

    :param list_of_lists: A list of lists
    :return result: A string that has been flattened from a list of lists
    """
    result = list()
    for i in list_of_lists:
        if isinstance(i, list):
            result.extend(flatten(i))
        else:
            result.append(str(i))
    return result


def expand_pandas(string: List, df: pd.DataFrame, allow_missing=False) -> List:
    """Expand string following columns in the dataframe"""
    return set(
        flatten(
            [
                expand(string, **row._asdict(), allow_missing=allow_missing)
                for row in df.itertuples(False)
            ]
        )
    )


def get_tmp(large: bool = False) -> str:
    """Returns path to temporary folder

    :param large: large temp folder (NFS), instead of local (e.g. /tmp/)
    :return path: a string with path to temp folder
    """
    import tempfile

    if large:
        path = Path("temp") / "large_temp"
    else:
        if "PBS_JOBID" in os.environ:
            path = Path("/scratch/$PBS_JOBID")
        elif "SLURM_JOB_ID" in os.environ and "SLURM_TMPDIR" in os.environ:
            path = Path("$SLURM_TMPDIR")
        else:
            path = tempfile.tempdir

    return str(path) + "/"


### Config


def _item_or_sample(row, item):
    i = getattr(row, item, None)
    if pd.isnull(i):
        return getattr(row, "sample")
    return i


def is_activated(xpath):
    c = config
    for entry in xpath.split("/"):
        c = c.get(entry, {})
    return bool(c.get("activate", False))


### Samples


def get_rule_stats(rule_name):
    r = re.compile("^stats/")
    return set(filter(r.match, getattr(rules, rule_name).output))


### Units


def get_read_type_trim(read_type_map):
    if read_type_map == "pe":
        return ["R1", "R2"]
    elif read_type_map == "se":
        return ["R"]
    else:
        return read_type_map
