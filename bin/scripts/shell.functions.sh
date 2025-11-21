#!/bin/bash
# Create a function for waiting job to finish (avoid to use all memory of the group in the server)
# made by Giovanni Ladeira, July 12, 2024
# Thank you Giovanni
function wait_for_job_completion() {
    # Crete a local variable (job_id) to store the first argument passed to the function
    local job_id=$1
    local time_min=$2
    echo "          Waiting for job $job_id to complete..."
    # While loop only stops when break or exit commands are executed
    while true; do
        # Check if the job is still in the queue
        job_status=$(squeue -j $job_id -h -o %T 2>/dev/null)
        if [[ -z $job_status ]]; then
            echo ""
            break
        elif [[ $job_status == "FAILED" ]]; then
            echo "          Job $job_id failed."
            exit 1
        else
            # Job is still running or pending
            echo "          Job $job_id is $job_status. Checking again in ${time_min} minutes..."
            sleep ${time_min}m
        fi
    done
}

#
# wait_for_job_completion_less_output.sh
#
# Usage:
#   wait_for_job_completion_less_output <job_id> [<check_interval_min>] [<message_interval_min>]
#
#   <job_id>             Slurm job ID to watch (required)
#   <check_interval_min> How many minutes between each squeue poll (default: 1)
#   <message_interval_min>
#                        After this many minutes of still‐running, print a status line (default: 30)
#
# Example:
#   wait_for_job_completion_less_output 123456
#     → polls every 1 minute; prints status every 30 minutes if job is still running

function wait_for_job_completion_less_output() {
    # --- Argument parsing & sanity checks ---

    # Must have at least 1 argument (job_id)
    if [[ $# -lt 1 ]]; then
        echo "Usage: wait_for_job_completion_less_output <job_id> [<check_interval_min>] [<message_interval_min>]"
        return 2
    fi

    local job_id="$1"
    # Check that job_id is numeric
    if ! [[ $job_id =~ ^[0-9]+$ ]]; then
        echo "ERROR: <job_id> must be a numeric Slurm job ID."
        return 2
    fi

    # Default check interval = 1 minute, message interval = 30 minutes
    local check_min="${2:-1}"
    local message_min="${3:-30}"

    # Validate that check_min is a positive integer
    if ! [[ $check_min =~ ^[0-9]+$ ]] || (( check_min <= 0 )); then
        echo "ERROR: <check_interval_min> must be a positive integer (minutes)."
        return 2
    fi

    # Validate that message_min is a positive integer
    if ! [[ $message_min =~ ^[0-9]+$ ]] || (( message_min <= 0 )); then
        echo "ERROR: <message_interval_min> must be a positive integer (minutes)."
        return 2
    fi

    # --- Initial sanity‐check echo ---
    echo "Watching job ID: $job_id, polling every $check_min minute(s), status messages every $message_min minute(s)."

    # Convert minutes to seconds
    local check_sec=$(( check_min * 60 ))
    local message_sec=$(( message_min * 60 ))

    # Track how many seconds have elapsed since last printed message
    local elapsed_since_message=0

    # --- Main polling loop ---
    while true; do
        # Query squeue for this job's status (e.g. RUNNING, PENDING, COMPLETED, FAILED, etc.)
        job_status=$(squeue -j "$job_id" -h -o %T 2>/dev/null)

        # If squeue returns empty → job is no longer in the queue → assume it finished
        if [[ -z "$job_status" ]]; then
            echo "Job $job_id has completed."
            return 0
        fi

        # If Slurm says FAILED, bail out immediately
        if [[ "$job_status" == "FAILED" ]]; then
            echo "Job $job_id failed (status=FAILED)."
            return 1
        fi

        # Otherwise → job is still pending/running
        # Only print once every message_interval
        if (( elapsed_since_message >= message_sec )); then
            echo "Job $job_id is $job_status. Checking again in ${check_min} minute(s)..."
            elapsed_since_message=0
        fi

        # Sleep for check_sec, then accumulate elapsed time
        sleep "$check_sec"
        (( elapsed_since_message += check_sec ))
    done
}
# Function to calculate heritability
calculate_heritability() {
    local var_g=$1
    local var_e=$2
    echo "scale=4; $var_g / ($var_g + $var_e)" | bc
}


function wait_job_check_progress_GWAS() {
    # Input arguments
    local job_file=$1
    local time_min=$2
    local start_time_gwas=$3
    local result_files=$4

    # Main loop to track job statuses
    while true; do
        # Read job IDs from the file
        local job_ids=($(cat "$job_file"))
        echo "Checking jobs: ${job_ids[@]}"
        
        # Track incomplete jobs
        local remaining_jobs=()
        
        for job_id in "${job_ids[@]}"; do
            # Check if the job is still in the queue
            job_status=$(squeue -j "$job_id" -h -o %T 2>/dev/null)
            if [[ -n $job_status ]]; then
                if [[ $job_status == "FAILED" ]]; then
                    echo "          Job $job_id failed."
                    exit 1
                else
                    remaining_jobs+=("$job_id")
                fi
            fi
        done

        # If there are no remaining jobs, break the loop
        if [[ ${#remaining_jobs[@]} -eq 0 ]]; then
            echo "All jobs completed."
            break
        fi

        # Update job file with remaining jobs
        echo "${remaining_jobs[@]}" > "$job_file"
        echo "Jobs still running: ${remaining_jobs[@]}"
        
        # Capture the intermediate time after job check
        check_time_gwas=$(date +%s)
        
        # Calculate time difference from the start
        elapsed_seconds=$((check_time_gwas - start_time_gwas))
        
        # Convert seconds to days, hours, minutes, and seconds
        days=$((elapsed_seconds / 86400))
        hours=$(( (elapsed_seconds % 86400) / 3600 ))
        minutes=$(( (elapsed_seconds % 3600) / 60 ))
        seconds=$((elapsed_seconds % 60))
        
        # Calculate the total number of lines in all files within the directory
        num_lines=$(wc -l "${result_files}"/* 2>/dev/null | tail -1 | awk '{print $1}')
        
        # Count the number of files in the directory
        num_files=$(ls -1 "${result_files}" 2>/dev/null | wc -l)
        
        # Calculate the number of SNPs by subtracting the number of files from the total number of lines
        num_snps=$((num_lines - num_files))
        
        # Display the number of SNPs processed so far and the elapsed time
        echo "$num_snps SNPs done in ${days} days, ${hours} hours, ${minutes} minutes, ${seconds} seconds"
        
        # Wait before checking again
        echo "Checking again in ${time_min} minutes..."
        sleep "${time_min}m"
    done
}



function wait_for_jobs_array() {
    # --- 1) PARSE & VALIDATE INPUTS ---
    if [[ $# -lt 1 ]]; then
        echo "Usage: wait_for_jobs_completion_less_output <array_name> [<check_interval_min>] [<message_interval_min>]"
        return 1
    fi

    local job_array_name="$1"
    local check_min="${2:-1}"
    local message_min="${3:-30}"

    # Verify that check_min and message_min are positive integers
    if ! [[ "$check_min" =~ ^[0-9]+$ ]] || (( check_min <= 0 )); then
        echo "ERROR: <check_interval_min> must be a positive integer."
        return 1
    fi
    if ! [[ "$message_min" =~ ^[0-9]+$ ]] || (( message_min <= 0 )); then
        echo "ERROR: <message_interval_min> must be a positive integer."
        return 1
    fi

    # Confirm that the named variable exists and is an indexed array
    local declare_output
    if ! declare_output="$(declare -p "$job_array_name" 2>/dev/null)"; then
        echo "ERROR: Array '$job_array_name' does not exist."
        return 1
    fi
    if [[ ! "$declare_output" =~ "declare -a" ]]; then
        echo "ERROR: '$job_array_name' is not an indexed array."
        return 1
    fi

    # Use a nameref to access the actual array of job IDs
    local -n job_ids_ref="$job_array_name"

    # Ensure the array is non-empty
    if (( ${#job_ids_ref[@]} == 0 )); then
        echo "ERROR: Array '$job_array_name' is empty."
        return 1
    fi

    # --- 2) INITIAL OUTPUT (“SANITY CHECK”) ---
    echo "Polling every ${check_min} minute(s), status messages every ${message_min} minute(s)."
    echo "Watching Jobs:"
    {
        local count=0
        for job in "${job_ids_ref[@]}"; do
            printf "%s " "$job"
            (( count++ ))
            if (( count % 10 == 0 )); then
                echo
            fi
        done
        # If we didn’t end on a multiple of 10, print a final newline
        if (( count % 10 != 0 )); then
            echo
        fi
    }

    # --- 3) PREPARE FOR POLLING LOOP ---
    # Copy all job IDs into a “remaining” array; as each job finishes, we’ll remove it.
    local remaining_jobs=( "${job_ids_ref[@]}" )
    # How many seconds between each squeue check
    local check_sec=$(( check_min * 60 ))
    # How many seconds between each “status message” window
    local message_sec=$(( message_min * 60 ))
    # Track elapsed seconds since the last printed “status message”
    local elapsed_since_message=0

    # --- 4) POLLING LOOP ---
    while true; do
        # 4.1) Sleep one check‐interval
        sleep "$check_sec"
        (( elapsed_since_message += check_sec ))

        # 4.2) Determine which of the “remaining_jobs” have now finished
        local -a newly_finished=()
        for jid in "${remaining_jobs[@]}"; do
            # Query Slurm: if squeue -j <jid> yields empty, it’s no longer in the queue
            if [[ -z "$(squeue -j "$jid" -h -o %T 2>/dev/null)" ]]; then
                newly_finished+=( "$jid" )
            fi
        done

        # 4.3) Remove any newly finished jobs from remaining_jobs
        if (( ${#newly_finished[@]} > 0 )); then
            local -a updated_remaining=()
            for jid in "${remaining_jobs[@]}"; do
                local is_new=false
                for nf in "${newly_finished[@]}"; do
                    if [[ "$jid" == "$nf" ]]; then
                        is_new=true
                        break
                    fi
                done
                if ! $is_new; then
                    updated_remaining+=( "$jid" )
                fi
            done
            remaining_jobs=( "${updated_remaining[@]}" )
        fi

        # 4.4) If enough time has passed since the last status-message window,
        #       and there are jobs that newly finished, print them now
        if (( elapsed_since_message >= message_sec )) && (( ${#newly_finished[@]} > 0 )); then
            echo "Status message at ${message_min}min"
            echo "New Jobs Finished:"
            {
                local cnt=0
                for finished in "${newly_finished[@]}"; do
                    printf "%s " "$finished"
                    (( cnt++ ))
                    if (( cnt % 10 == 0 )); then
                        echo
                    fi
                done
                if (( cnt % 10 != 0 )); then
                    echo
                fi
            }
            echo
            elapsed_since_message=0
        fi

        # 4.5) If no jobs remain, exit the loop
        if (( ${#remaining_jobs[@]} == 0 )); then
            break
        fi
        # Otherwise, loop again (sleep → check → maybe print → etc.)
    done

    # --- 5) ALL DONE ---
    echo "All Jobs Finished"
    return 0
}
