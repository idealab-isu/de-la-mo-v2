# Copyright 2016-2018 Iowa State University Research Foundation, Inc.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


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

