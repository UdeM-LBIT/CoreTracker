import sys


class Output(object):
    """A simple Class to output results, either in a file or to stdout"""

    def __init__(self, file=None):
        self.file = file

        if(self.file):
            out = open(file, 'w')
            self.out = out
        else:
            self.out = sys.stdout

    def write(self, line):
        """ Write output """
        self.out.write(line)

    def close(self):
        """ Close Stream """
        if self.out is not sys.stdout:
            self.out.close()

    @staticmethod
    def error(message, type="Error"):
        sys.stderr.write("%s: %s\n" % (type, message))
