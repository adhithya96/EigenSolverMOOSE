#* This file is part of the MOOSE framework
#* https://www.mooseframework.org
#*
#* All rights reserved, see COPYRIGHT for full restrictions
#* https://github.com/idaholab/moose/blob/master/COPYRIGHT
#*
#* Licensed under LGPL 2.1, please see LICENSE for details
#* https://www.gnu.org/licenses/lgpl-2.1.html

from RunApp import RunApp
import os, signal,time, platform, subprocess, copy
from TestHarness import util
from tempfile import SpooledTemporaryFile


# Classes that derive from this class are expected to write
# output files. The Tester::getOutputFiles() method should
# be implemented for all derived classes.
class SignalTester(RunApp):
    @staticmethod
    def validParams():
        params = RunApp.validParams()
        params.addParam('signal', "SIGUSR1", "The signal to send to the app. Defaults to SIGUSR1")
        params.addParam('sleep_time', 1, "The amount of time the tester should wait before sending a signal.")
        return params

    def __init__(self, name, params):
        RunApp.__init__(self, name, params)

    def send_signal(self,pid):
        """Function used to send a signal to the program automatically for testing purposes."""

        #Create a while loop that checks if the stdout buffer has any data in it, and then sends the signal once
        #it knows that the moose_test binary is actually doing something.

        #process.poll() will return None if the process is running in the OS.
        #This acts as a safety precaution against an infinite loop -- this will always close.
        while(self.process.poll() == None):

            #first, make a true duplicate of the stdout file so we don't mess with the seek on the actual file
            out_dupe = copy.copy(self.outfile)
            #go to the beginning of the file and see if its actually started running the binary
            out_dupe.seek(0)
            output = out_dupe.read()

            #if the output is blank, the moose_test binary hasn't actually started doing anything yet.
            #if so, sleep briefly and check again.
            if not output:
                time.sleep(0.05)
                continue

            #if the output isn't blank, then we actually sleep for the time specified in sleep_time
            #then we finally send the SIGUSR1 and exit the loop
            time.sleep(self.specs['sleep_time'])
            os.kill(pid,signal.SIGUSR1)
            break

    def runCommand(self, timer, options):
        """
        Helper method for running external (sub)processes as part of the tester's execution.  This
        uses the tester's getCommand and getTestDir methods to run a subprocess.  The timer must
        be the same timer passed to the run method.  Results from running the subprocess is stored
        in the tester's output and exit_code fields.
        """

        exit_code = super().spawnSubprocessFromOptions(timer, options)
        if exit_code: # Something went wrong
            return

        self.send_signal(self.process.pid)
        super().finishAndCleanupSubprocess(timer)
