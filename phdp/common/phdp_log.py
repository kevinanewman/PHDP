"""

**Functions to manage batch and session log files.**

----

**CODE**

"""

print('importing %s' % __file__)

from common import phdp_globals  # import global variables


def init_logfile():
    """
    Create a session logfile.

    """
    import time, datetime
    import common.file_io as file_io
    from phdp import code_version

    file_io.validate_folder(phdp_globals.options.output_folder)

    phdp_globals.options.logfilename = '%s%s.txt' % (
        phdp_globals.options.output_folder + phdp_globals.options.logfile_prefix,
        phdp_globals.options.session_unique_name)

    with open(phdp_globals.options.logfilename, 'w') as log:
        log.write('PHDP %s session %s started at %s %s\n\n' % (
            code_version, phdp_globals.options.session_name, datetime.date.today(), time.strftime('%H:%M:%S')))


def end_logfile(message):
    """
    End logfile with closing message, record elapsed time.

    Args:
        message (str): message string to write

    """
    import time
    phdp_globals.options.end_time = time.time()
    elapsed_time = (phdp_globals.options.end_time - phdp_globals.options.start_time)
    import datetime
    logwrite('\nSession ended at %s %s' % (datetime.date.today(), time.strftime('%H:%M:%S')))
    logwrite('Session elapsed time %.2f seconds\n' % elapsed_time)
    logwrite(message, terminator='')


def logwrite(message, echo_console=True, terminator='\n'):
    """
    Write message to logfile.

    Args:
        message (str or [strs]): message string or list of strings to write
        echo_console (bool): write message to console if True
        terminator (str): end of message terminator, default is newline (``\\n``)

    """
    with open(phdp_globals.options.logfilename, 'a') as log:
        if type(message) is list:
            for m in message:
                log.write(m + terminator)
        else:
            log.write(message + terminator)
        if phdp_globals.options.verbose or echo_console:
            print(message)
