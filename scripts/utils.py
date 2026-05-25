"""Shared utilities for genDiv analysis scripts."""

import pandas as pd
import numpy as np

# Shared constants
BOOK_SIZE_ORDER = ['HIGH', 'MEDIUM', 'LOW']
BOOK_SIZE_COLORS = {'HIGH': '#E74C3C', 'MEDIUM': '#F39C12', 'LOW': '#2ECC71'}


def format_stats(group_df, columns):
    """Return a dict of 'mean +/- SD' strings for the given columns."""
    results = {}
    for col in columns:
        m = group_df[col].mean()
        s = group_df[col].std()
        results[col] = f"{m:.4f} +/- {s:.4f}"
    return results
