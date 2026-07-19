#!/bin/bash
SCRIPT=$1
START=$2
END=$3
CHUNK=10

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

    while squeue -j ${JOBID} -h 2>/dev/null | grep -q .; do
        sleep 15
    done
    echo "Range ${current}-${stop} complete"

    current=$((stop + 1))
done

echo "All requested ranges submitted and completed"
