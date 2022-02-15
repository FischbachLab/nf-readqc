#!/usr/bin/env python3

import logging
import sys
from subprocess import run


def concatenate(fastq_list, output_name):
    all_fastq_joined = " ".join(fastq_list)
    cmd = f"zcat {all_fastq_joined} | gzip > {output_name}"
    logging.info(f"Running: {cmd}")
    run(cmd, shell=True)
    return


def all_same_direction(fastq_list, fastq_type):
    same_direction = all(
        [True if fastq_type in fastq else False for fastq in fastq_list]
    )
    if same_direction:
        return True
    else:
        logging.error(
            f"Cannot concatenate fastq files of different orientation, expected all to be of type {fastq_type}, got {fastq_list}"
        )
        return False


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.DEBUG, format="%(asctime)s %(levelname)s %(message)s"
    )
    prefix = sys.argv[1]
    single_end = sys.argv[2]  # true/false
    list_of_fastqs = sorted(sys.argv[3:])

    se_bool = True if single_end == "true" else False

    cat_fwd = f"{prefix}.merged.R1.fastq.gz"
    cat_rev = f"{prefix}.merged.R2.fastq.gz"

    logging.info(f"Is SE?: {se_bool}")
    if se_bool:
        assert all_same_direction(list_of_fastqs, "R1")
        concatenate(list_of_fastqs, cat_fwd)
    else:
        # if paired end is true, a sorted list will contain
        # fwd and rev in interleaved order (fwd_1, rev_1, fwd_2, rev_2, ...)

        # extract all even items from the list (0, 2, 4 ...)
        all_fwd = list_of_fastqs[::2]
        assert all_same_direction(all_fwd, "R1")
        logging.info(f"Concatenating all FWD: {all_fwd}")
        concatenate(all_fwd, cat_fwd)

        # extract all even items from the list (1, 3, 5 ...)
        all_rev = list_of_fastqs[1::2]
        assert all_same_direction(all_rev, "R2")
        logging.info(f"Concatenating all REV: {all_rev}")
        concatenate(all_rev, cat_rev)
