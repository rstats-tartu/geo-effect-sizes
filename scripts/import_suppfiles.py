import os
import re
import gzip
import tarfile
import io
import collections
import argparse
import pandas as pd
import numpy as np
from scipy.stats import binom
from pandas.api.types import is_string_dtype
import numbers
import warnings


class FormatError(Exception):
    pass


xls = re.compile("xls")
drop = "series_matrix\\.txt\\.gz$|filelist\\.txt$|readme|\\.bam(\\.tdf|$)|\\.bai(\\.gz|$)|\\.sam(\\.gz|$)|\\.csfasta|\\.fa(sta)?(\\.gz|\\.bz2|\\.txt\\.gz|$)|\\.f(a|n)a(\\.gz|$)|\\.wig|\\.big[Ww]ig$|\\.bw(\\.|$)|\\.bed([Gg]raph)?(\\.tdf|\\.gz|\\.bz2|\\.txt\\.gz|$)|(broad_)?lincs|\\.tdf$|\\.hic$|\\.rds(\\.gz|$)|\\.tar\\.gz$|\\.mtx(\\.gz$|$)|dge\\.txt\\.gz$|umis?\\.txt\\.gz$"
drop = re.compile(drop)
pv_str = "p[^a-zA-Z]{0,4}val"
pv = re.compile(pv_str)
adj = re.compile("adj|fdr|corr|thresh|q[^a-zA-Z]{0,4}val")
fc = re.compile("(l(og)?([0-9])?((\\()?f(old)?(_)?c|ratio)|^fc$)")
ws = re.compile(" ")
mtabs = re.compile("\\w+\\t{2,}\\w+")
tab = re.compile("\\t")
fields = ["es_var", "es_val", "pv_var", "pv_val", "adjpv_var", "adjpv_val", "note"]
PValSum = collections.namedtuple("PValSum", fields, defaults=[np.nan] * len(fields))
peak = re.compile("(narrow|broad|gapped)peak")


class ImportSuppfiles(object):
    def __init__(self):
        self.out = {}

    def from_flat(self, input, tar=None):
        if drop.search(input.name.lower() if tar else input.lower()):
            key = os.path.basename(input.name if tar else input)
            return self.out.update(note(key, "not imported"))
        else:
            out = {}
            try:
                if xls.search(input.name if tar else input):
                    try:
                        out.update(self.read_excel(input, tar=tar))
                    except ValueError as e:
                        out.update(self.read_csv(input, tar=tar))
                else:
                    d = self.read_csv(input, tar=tar)
                    is_empty = [v.empty for v in d.values()][0]
                    if is_empty:
                        raise Exception("empty table")
                    else:
                        peakfile = peak.search(
                            input.name.lower() if tar else input.lower()
                        )
                        if peakfile:
                            key = os.path.basename(input.name if tar else input)
                            d[key].loc[-1] = d[key].columns
                            d[key] = d[key].sort_index().reset_index(drop=True)
                            d[key].columns = eval(peakfile.group(0))
                        out.update(d)
            except Exception as e:
                key = os.path.basename(input.name if tar else input)
                peakfile = peak.search(input.name.lower() if tar else input.lower())
                if peakfile:
                    e = f"Misspecified '{peakfile.group(0)}' file; {e}"
                out.update(note(key, e))
            return self.out.update(out)

    def from_tar(self, input):
        with tarfile.open(input, "r:*") as tar:
            for member in tar:
                if member.isfile():
                    self.from_flat(member, tar)

    def find_header(self, df, n=20):
        head = df.head(n)
        matches = [
            i[0]
            for i in [
                [i for i, x in enumerate(head[col].str.contains(pv_str, na=False)) if x]
                for col in head
            ]
            if i
        ]
        idx = min(matches) + 1 if matches else 0
        if idx == 0:
            for index, row in head.iterrows():
                if all([isinstance(i, str) for i in row if i is not np.nan]):
                    idx = index + 1
                    break
        return idx

    def csv_helper(self, input, input_name, csv, verbose=0):
        # Get comments and set rows to skip
        r = pd.read_csv(csv, sep=None, engine="python", iterator=True, nrows=1000)
        comment = None
        sep = r._engine.data.dialect.delimiter
        columns = r._engine.columns
        if isinstance(input, (tarfile.ExFileObject)):
            with csv as h:
                first_line = h.readline()
        elif input_name.endswith("gz") or isinstance(input, (gzip.GzipFile)):
            with gzip.open(input) as h:
                first_line = h.readline().decode("utf-8").rstrip()
        else:
            with open(input, "r") as h:
                first_line = h.readline().rstrip()
        more_tabs_than_sep = len(tab.findall(first_line)) > len(
            re.findall(sep, first_line)
        )
        if re.search("^#", first_line) or more_tabs_than_sep:
            comment = "#"
            # Get delimiter
            r = pd.read_csv(
                csv, sep=None, engine="python", iterator=True, skiprows=20, nrows=1000
            )
            sep = r._engine.data.dialect.delimiter
            columns = r._engine.columns
        if ws.search(sep):
            sep = "\\s+"
        if mtabs.search(first_line):
            sep = "\\t+"
        # Import file
        if isinstance(input, (tarfile.ExFileObject)) and input_name.endswith("gz"):
            with gzip.open(input) as h:
                df = pd.read_csv(h, sep=sep, comment=comment, encoding="unicode_escape")
        else:
            df = pd.read_csv(input, sep=sep, comment=comment, encoding="unicode_escape")
        # Check and fix column names
        # Case of extra level of delimiters in column names
        if len(df.columns) > len(columns):
            df = pd.read_csv(
                input,
                header=None,
                skiprows=[0],
                sep=sep,
                comment=comment,
                encoding="unicode_escape",
            ).drop([0])
            df.columns = columns
        unnamed = ["Unnamed" in i for i in df.columns]
        # Case of empty rows before header
        if all(unnamed):
            idx = self.find_header(df)
            if idx > 0:
                df = pd.read_csv(
                    input,
                    sep=sep,
                    comment=comment,
                    skiprows=idx,
                    encoding="unicode_escape",
                )
        # Case of anonymous row names
        if unnamed[-1] & sum(unnamed) == 1:
            if any([pv.search(i) for i in df.columns]):
                df.columns = [df.columns[-1]] + list(df.columns[:-1])
        if verbose > 1:
            print("df after import:\n", df)
        return {os.path.basename(input_name): df}

    def excel_helper(self, input, input_name, verbose=0):
        tabs = {}
        if input_name.endswith("gz") or isinstance(input, (gzip.GzipFile)):
            excel_file = gzip.open(input)
        else:
            excel_file = input
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            wb = pd.ExcelFile(excel_file)
            if len(wb.sheet_names) == 0:
                (m,) = [i.message for i in w][0].args
                raise FormatError(
                    f"The data source could not be successfully parsed with warning: '{m}'",
                )
        sheets = wb.sheet_names
        sheets = [i for i in sheets if "README" not in i]
        for sheet in sheets:
            df = wb.parse(sheet, comment="#")
            if df.empty:
                df = wb.parse(sheet)
            if verbose > 1:
                print("df after import:\n", df)
            if not df.empty:
                pu = sum(["Unnamed" in i for i in list(df.columns)]) / len(df.columns)
                if pu >= 2 / 3:
                    idx = self.find_header(df)
                    if idx > 0:
                        df = wb.parse(sheet, skiprows=idx)
                tabs.update({os.path.basename(input_name) + "-sheet-" + sheet: df})
        return tabs

    def read_csv(self, input, tar=None):
        if isinstance(input, (tarfile.TarInfo)):
            input_name = os.path.basename(input.name)
            with tar.extractfile(input) as h:
                if input_name.endswith("gz"):
                    with gzip.open(h) as gz:
                        csv = io.StringIO(gz.read().decode("unicode_escape"))
                else:
                    csv = io.StringIO(h.read().decode("unicode_escape"))
            with tar.extractfile(input) as h:
                out = self.csv_helper(h, input_name, csv)
        else:
            input_name = input
            csv = input
            out = self.csv_helper(input, input_name, csv)
        return out

    def read_excel(self, input, tar=None):
        if isinstance(input, (tarfile.TarInfo)):
            input_name = os.path.basename(input.name)
            with tar.extractfile(input) as h:
                out = self.excel_helper(h, input_name)
        else:
            input_name = input
            out = self.excel_helper(input, input_name)
        return out


def get_pvalues(i):
    return bool(pv.search(i.lower()))


def get_adj_pvalues(i):
    return bool(adj.search(i.lower()))


def get_fc(i):
    il = i.lower()
    return bool(fc.search(il) and not re.search("lfcse", il))


def fix_column_dtype(df):
    for col in df.columns:
        s = df[col]
        if is_string_dtype(s):
            if "," in s[:5].astype(str).str.cat(sep=" "):
                df[col] = s.apply(lambda x: str(x).replace(",", "."))
            df[col] = pd.to_numeric(s, errors="coerce")
    return df


def split_extra_delim(df, cols, delim, func):
    # Check if there is ANOTHER(fucking!!#?) level of ":" delimiters in column names
    for index, col in enumerate(cols):
        col_count = len(re.findall(delim, col))
        obs_count = len(re.findall(delim, str(df.iloc[0, index])))
        if obs_count == 0:
            pass
        elif col_count == obs_count:
            new_cols = col.split(delim)
            split_pval_col = [i for i in new_cols if func(i)]
            cols_split = df.iloc[:, index].str.split(delim, expand=True)
            try:
                cols_split.columns = new_cols
                df[split_pval_col] = cols_split[split_pval_col]
                df.drop(col, axis=1, inplace=True)
            except ValueError:
                pass
    return df


def check_extra_delim(df, delim, func):
    split_cols = [i for i in df.columns if delim in i]
    if split_cols:
        df = split_extra_delim(df, cols=split_cols, delim=delim, func=func)
    return df


def parse_table(df):

    # Drop columns with numeric column names
    df = df.filter(regex="^\\D")
    # Drop columns with NaN column names
    df = df.loc[:, df.columns.notnull()]
    df.columns = map(str.lower, df.columns)
    pval_cols = [i for i in df.columns if get_pvalues(i)]
    adj_cols = [i for i in df.columns if get_adj_pvalues(i)]
    fc_cols = [i for i in df.columns if get_fc(i)]
    # some outputs have both fold change (fc) and log fold change columns, let's keep only logFC
    fold_changes = [i for i in fc_cols if re.search("^fc|(?<!l)fc$", i)]
    if fold_changes and len(fold_changes) < len(fc_cols):
        fc_cols = [i for i in fc_cols if i not in fold_changes]
    pvalues = df[pval_cols].copy()
    adjpval = df[adj_cols].copy()
    fcval = df[fc_cols].copy()
    # Check if there is ANOTHER(!!#?) level of ":" delimiters in p value column(s)
    pvalues = check_extra_delim(pvalues, delim=":", func=get_pvalues)
    adjpval = check_extra_delim(adjpval, delim=":", func=get_adj_pvalues)
    fcval = check_extra_delim(fcval, delim=":", func=get_fc)
    # fix data type
    pvalues_check = fix_column_dtype(pvalues)
    adjpval_check = fix_column_dtype(adjpval)
    fcval_check = fix_column_dtype(fcval)
    merged = (
        fcval_check.melt(var_name="es_var", value_name="es_val")
        .join(pvalues_check.melt(var_name="pv_var", value_name="pv_val"))
        .join(adjpval_check.melt(var_name="adjpv_var", value_name="adjpv_val"))
    )
    return merged


def note(filename, message):
    return {
        filename: pd.DataFrame(PValSum(note=str(message).rstrip())._asdict(), index=[0])
    }


def parse_key(k, filename):
    if k == filename:
        key = k
    else:
        replacement = r"\2" if filename in k else r"\2 from \1"
        key = re.sub(r"^(.*)-(sheet-.*)", replacement, k) + " from " + filename
    return key


def parse_suppfiles(
    file=None,
    list=None,
    blacklist=None,
    verbose=False,
    **kwargs,
):
    def single_true(iterable):
        i = iter(iterable)
        return any(i) and not any(i)

    assert single_true(
        [file, list]
    ), "Provide either supplementary file name (file) or file with supplementary file names (list)."

    if file:
        input = file
    elif list:
        input = []
        with list as f:
            for line in f:
                input.append(line.rstrip())

    blacklist = []
    if blacklist:
        with blacklist as f:
            for line in f:
                file = os.path.basename(line.rstrip())
                blacklist.append(file)

    # Keep only inputs that exist
    input = [i for i in input if os.path.isfile(i)]

    # Drop files in blacklist
    if blacklist:
        input = [i for i in input if os.path.basename(i) not in blacklist]

    assert len(input) > 0, "No files to process."

    if verbose:
        print("Working on ", input)

    res = {}
    for path in input:
        filename = os.path.basename(path)
        frames = ImportSuppfiles()
        if tarfile.is_tarfile(path):
            frames.from_tar(path)
        else:
            frames.from_flat(path)
        frames = {
            k: v
            if all(i in fields for i in v.columns)
            else parse_table(v)
            if any(
                [get_pvalues(i) for i in v.columns if not isinstance(i, numbers.Number)]
            )
            else pd.DataFrame(PValSum(note="no pvalues")._asdict(), index=[0])
            for k, v in frames.out.items()
        }
        for k, v in frames.items():
            res.update({parse_key(k, filename): v})

    result = pd.concat(
        [df for df in res.values()],
        keys=[k for k in res.keys()],
        names=["id"],
        sort=False,
    )

    return result


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--file", nargs="+", help="path to input file to be parsed")
    group.add_argument(
        "--list",
        metavar="FILE",
        type=argparse.FileType("r"),
        help="file with paths to input files, one per line",
    )
    parser.add_argument("--out", metavar="FILE", help="output file", required=True)
    parser.add_argument(
        "--verbose", "-v", help="increase output verbosity", action="store_true"
    )
    parser.add_argument(
        "--blacklist",
        metavar="FILE",
        type=argparse.FileType("r"),
        help="file with filenames to skip importing, one per line",
    )
    args = parser.parse_args()

    df = parse_suppfiles(**args.__dict__)

    with open(args.out, "w") as f:
        df.reset_index(level="id").to_csv(f, index=False)
