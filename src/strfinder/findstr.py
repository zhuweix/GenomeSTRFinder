import re
from typing import Annotated
from pathlib import Path
import tomllib
import pysam
import typer
import pandas as pd
import numpy as np

def count_str_repeats(seq: str, motif_pat: re.Pattern, 
                      motif_size: int, max_dist: int) -> list:
    def count_motif(match):
        start = match.start()
        end = match.end()
        count = (end - start) // motif_size
        return (start, end, count)

    # combine adjacent motifs
    match_list = [count_motif(match) for match in re.finditer(motif_pat, seq)]
    if len(match_list) == 0:
        return []
    cur_motif = match_list[0]
    _, ce, cc = cur_motif
    cc = [cc, ]
    if len(match_list) == 1:
        return [str(cc[0])]
    result = []
    for (s, e, c) in match_list[1:]:
        if s - ce <= max_dist:
            ce = e
            cc.append(c)
        else:
            result.append(';'.join(map(str, cc)))
            ce = e
            cc = [c]
    result.append(';'.join(map(str, cc)))
    return result

def build_pattern(motif: str, min_count: int) -> tuple:
    # Escape any special regex characters in the motif
    escaped_motif = re.escape(motif)
    # Build the pattern string
    pattern = f"(?:{escaped_motif}){{{min_count},}}"
    # Compile and return the pattern
    return re.compile(pattern)

def process_str_fn(str_fn: str, min_count: int) -> dict:
    str_pd = pd.read_csv(str_fn, index_col=False)
    std_dict = {}
    for chrom, tmp_pd in str_pd.groupby('chrom', sort=False):
        std_dict[chrom] = [(name, motif, build_pattern(motif, min_count)) 
                           for (name, motif) in zip(tmp_pd['name'], tmp_pd['seq'])]
    return std_dict

def transform_ref(chrom_fn: str):
    chrom_pd = pd.read_csv(chrom_fn, index_col=False)
    chrom_dict = dict(zip(chrom_pd['ref'], chrom_pd['chrom']))
    return chrom_dict


def count_str_bam(bam_fn: str, motif_fn: str, 
                  min_count: int, chrom_fn: str, max_dist: int) -> pd.DataFrame:
    chrom_dict = transform_ref(chrom_fn)
    motif_dict = process_str_fn(motif_fn, min_count)
    result_pd = []
    with pysam.AlignmentFile(bam_fn, 'rb') as inbam:
        for read in inbam:
            ref = read.reference_name
            if ref is None:
                continue
            chrom = chrom_dict.get(ref, None)
            if chrom is None:
                continue
            qseq = read.query_alignment_sequence
            if qseq is None:
                continue
            motif_list = motif_dict.get(chrom, None)
            if motif_list is None:
                continue
            for (name, motif, pat) in motif_list:
                count = count_str_repeats(qseq, pat, len(motif), max_dist)
                if len(count) > 0:
                    tmp_pd = pd.DataFrame({'STRName': name, 'Motif': motif, 'MotifCount': count})
                    result_pd.append(tmp_pd)
    result_pd = pd.concat(result_pd, ignore_index=True)
    return result_pd



def process_bam(
    bam_fn: str, config_fn: str, output_fn: str):
    with open(config_fn, 'rb') as filep:
        config_dict = tomllib.load(filep)
    chrom_fn = config_dict['Data']['chrom_file']
    str_fn = config_dict['Data']['str_file']
    min_count = config_dict['Parameter']['min_count']
    max_dist = config_dict['Parameter']['max_dist']
    result_pd = count_str_bam(
        bam_fn, str_fn, min_count,
        chrom_fn, max_dist)
    result_pd.to_csv(output_fn, index=False)

app = typer.Typer(help="Process BAM file and count STRs based on configuration.")
@app.command()
def cli(
    bam_fn: Annotated[
        Path,
        typer.Argument(
            exists=True,
            dir_okay=False,
            help="Path to input BAM file"
        )
    ],
    config_fn: Annotated[
        Path,
        typer.Argument(
            exists=True,
            dir_okay=False,
            help="Path to TOML configuration file"
        )
    ],
    output_fn: Annotated[
        Path,
        typer.Argument(
            dir_okay=False,
            help="Path to output CSV file"
        )
    ]):
    process_bam(bam_fn, config_fn, output_fn)

if __name__ == "__main__":
    app()