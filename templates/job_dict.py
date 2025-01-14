main_job = {
    'universe': "vanilla",
    'executable': "run.sh",
    'JobBatchName': "batchname", 
    '+singularityimage': '"/data9/Users/choij/Singularity/images/private-el9"',
    'RequestMemory': f'2000 MB', # 2 times of memory usage
    'RequestCpus': 1,
    'arguments': "$(Process)",
    'output': "job_$(Process).out",
    'error': "job_$(Process).err",
    'should_transfer_files': "YES",
    'when_to_transfer_output': "ON_EXIT_OR_EVICT",
    'concurrency_limits': "300",
    # periodic release for 3 times // Held reason is lack of memory // JobStatus is Hold // https://research.cs.wisc.edu/htcondor/manual/v8.5/12_Appendix_A.html
    'periodic_release': '(NumJobStarts < 5) && (HoldReasonCode == 34 || HoldReasonCode == 21) && (JobStatus == 5)' 
}

hadd_job = {
    'universe': "vanilla",
    'executable': "hadd.sh",
    'JobBatchName': "batchname",
    '+singularityimage': '"/data9/Users/choij/Singularity/images/private-el9"',
    'RequestCpus': 4,
    'RequestMemory': '8192 MB',
    'output': "hadd.out",
    'error': "hadd.err",
    'should_transfer_files': "YES",
    'when_to_transfer_output': "ON_EXIT_OR_EVICT",
    'periodic_release': '(NumJobStarts < 3) && (HoldReasonCode == 34 || HoldReasonCode == 21) && (JobStatus == 5)' # periodic release for 3 times // Held reason is lack of memeory // JobStatus is Hold // https://research.cs.wisc.edu/htcondor/manual/v8.5/12_Appendix_A.html
}

final_job = {
    'universe': "vanilla",
    'executable': "final.sh",
    'JobBatchName': "Summarize",
    'getenv': "True",
    'RequestMemory': '1024 MB',
    'output': "final.out",
    'error': "final.err",
    'should_transfer_files': "YES",
    'when_to_transfer_output': "ON_EXIT_OR_EVICT"
}