#! /usr/bin/env python3

import argparse
import pandas as pd


def fake_analysis(df, labels_col, label1, label2, deps_cols):
    col_names = df.columns.values.tolist()
    if labels_col not in col_names:
        raise ValueError(
            f"Cannot find labels column {labels_col} in the dataframe column names."
        )
    for dep_col in deps_cols:
        if dep_col not in col_names:
            raise ValueError(
                f"Cannot find dependent variable column {dep_col} in the dataframe column names."
            )
    return col_names


def get_command_line():
    ap = argparse.ArgumentParser()
    ap.add_argument("infile", help="CSV data file")
    return ap.parse_args()


if __name__ == "__main__":
    args = get_command_line()

    # Read the input data
    df = pd.read_csv(args.infile)

    labels_col = "Label"
    labels_col = "LABEL"
    deps_cols = ["DEP_001", "DEP_002", "DEP_003", "DEP_004", "DEP_005"]
    report = fake_analysis(df, labels_col, "Control", "Treated", deps_cols)
