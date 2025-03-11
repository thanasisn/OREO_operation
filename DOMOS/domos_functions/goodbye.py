

def goodbye(logfile, tic=tic, scriptname=SCRIPT_NAME, quiet=False):
    """
    Log script execution to a central file
    """
    out  = datetime.now().strftime("%F %T") + " "
    out += os.getlogin() + "@" + os.uname()[1] + " "
    out += scriptname    + " "
    out += str(round((datetime.now() - tic).total_seconds() / 60.0, 2)) + " mins"
    if not quiet:
        print('\n' + out + '\n')
    with open(logfile, 'a') as runlog:
        runlog.write(out + '\n')



