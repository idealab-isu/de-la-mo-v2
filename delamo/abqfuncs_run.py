
def RunJob(job,wait=False):
    # job is an abq.mdb.Job instance to run
    # if wait=True, the interface will sit there spinning
    # waiting for the job to complete
    # Otherwise job will run in the background
    job.submit()
    if wait:
        job.waitForCompletion()
        pass
    pass

