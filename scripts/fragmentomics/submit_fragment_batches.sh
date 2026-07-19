#!/bin/bash
SCRIPT=/scratch/alice/n/nhsas1/ptcl_ctdna_thesis/scripts/fragmentomics/run_fragment_extraction.slurm
CHUNK=10
START=$1
END=$2

current=${START}
while [ ${current} -le ${END} ]; do
    stop=$((current + CHUNK - 1))
    if [ ${stop} -gt ${END} ]; then
        stop=${END}
    fi

    echo "Submitting array ${current}-${stop}"
    JOBID=$(sbatch --array=${current}-${stop}%4 --parsable ${SCRIPT})

    if [ -z "${JOBID}" ]; then
        echo "ERROR: submission failed for range ${current}-${stop}, stopping"
        exit 1
    fi
    echo "Job ID: ${JOBID}"

    # poll until every task of this job has left the queue
    while squeue -j ${JOBID} -h 2>/dev/null | grep -q .; do
        sleep 15
    done
    echo "Range ${current}-${stop} complete"

    current=$((stop + 1))
done

echo "All requested ranges submitted and completed"
