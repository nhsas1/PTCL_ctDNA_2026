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

    # Leaving the queue is not the same as succeeding. A chunk that FAILED, was CANCELLED
    # or hit TIMEOUT disappears from squeue exactly as a completed one does, so without
    # this check the driver printed "complete" and submitted the next chunk regardless -
    # and finished by announcing that everything completed. A failed extraction chunk was
    # indistinguishable from a successful one in this script's own output.
    BAD=$(sacct -j "${JOBID}" --format=State -n -X 2>/dev/null \
          | tr -d ' ' | grep -v '^COMPLETED$' | sort -u | tr '\n' ' ')
    if [ -n "${BAD}" ]; then
        echo "ERROR: range ${current}-${stop} did not complete cleanly. States: ${BAD}"
        echo "Stopping so the failure is not carried forward silently."
        exit 1
    fi
    echo "Range ${current}-${stop} complete"

    current=$((stop + 1))
done

echo "All requested ranges submitted and completed"
