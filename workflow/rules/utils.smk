import pandas as pd
from typing import List


#################
### FUNCTIONS ###
#################

### General


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


def expand_pd(string: List, df: pd.DataFrame, allow_missing=False) -> List:
    return set(expand(string, zip, **df.to_dict("list"), allow_missing=allow_missing))

### Config


def is_activated(xpath):
    c = config
    for entry in xpath.split("/"):
        c = c.get(entry, {})
    return bool(c.get("activate", False))


### Samples


def get_rule_stats(rule_name):
    r = re.compile("^stats/")
    return set(filter(r.match, getattr(rules, rule_name).output))
