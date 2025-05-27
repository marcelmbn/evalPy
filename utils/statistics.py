"""
Python module for statistical functions.
This module contains functions for calculating the MAE and WTMAD-2.
"""

import numpy as np
import pandas as pd


def statistical_measures(subset_data: pd.DataFrame) -> dict[str, int | float]:
    """
    Calculate the statistics for a subset of the dataframe.
    """
    errors = subset_data["MethodValue"] - subset_data["ReferenceValue"]
    stats_dict = {
        "N": len(subset_data),
        "MeanAbsRef": np.mean(np.abs(subset_data["ReferenceValue"])),
        "MAE": np.mean(np.abs(errors)),
        "MSE": np.mean(errors),
        "STDDEV": np.std(errors, ddof=1),
        "RMSD": np.sqrt(np.mean(errors**2)),
        "MAX": np.max(errors),
        "MIN": np.min(errors),
        "ErrRange": np.max(errors) - np.min(errors),
    }
    return stats_dict
