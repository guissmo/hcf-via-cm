import atexit
import os
import time

# This file sets
# - DOT_SAGE
# - LOGFILEDIR
# - RECIPDIR
# - FGAGDIR
# - SHIMURADIR
# if they exist.

ip = get_ipython()

# Input history will have no limit.
LIMIT = 0

# Set DOT_SAGE variable in SAGE to the environment variable DOT_SAGE.
DOT_SAGE = "";
if 'DOT_SAGE' in os.environ.keys():
    DOT_SAGE = os.environ['DOT_SAGE'];

# Set LOGFILEDIR variable in SAGE to the environment variable LOGFILEDIR, if it is set.
if 'LOGFILEDIR' in os.environ.keys():
    LOGFILEDIR = os.environ['LOGFILEDIR'];
    sessionname = time.strftime('%Y%m%d-%H%M');
    F = time.strftime(LOGFILEDIR+'/%Y%m%d-%H%M.sagelog');
    D = time.strftime(LOGFILEDIR+'/%Y%m%d-%H%M.sagedump');
    histfile = os.path.join(os.path.expanduser(F));
    dumpfile = os.path.join(os.path.expanduser(D));
    print("Recording input history to %s" % histfile);
    print("Recording output history to %s" % dumpfile);
    load(DOT_SAGE+"/histsave.sage");
    atexit.register(save_history);
else:
    print("No log files.");

# Set RECIPDIR variable in SAGE to the environment variable RECIPDIR, if it is set.
if 'RECIPDIR' in os.environ.keys():
    RECIPDIR = os.environ['RECIPDIR'];
    load("recip.sage");
    print("Recip loaded from %s" % RECIPDIR);
else:
    RECIPURL = "https://bitbucket.org/mstreng/recip/raw/master/recip_online.sage";
    load(RECIPURL);
    print("RECIP loaded from %s" % RECIPURL);

# Set FGAGDIR variable in SAGE to the environment variable FGAGDIR, if it is set.
if 'FGAGDIR' in os.environ.keys():
    FGAGDIR = os.environ['FGAGDIR'];
    gp("FGAGDIR=\""+FGAGDIR+"\";");
    # print("FGAGDIR = %s" % FGAGDIR);
else:
    print("NO FGAGDIR.");

# Set SHIMURADIR variable in SAGE to the environment variable SHIMURADIR, if it is set.
if 'SHIMURADIR' in os.environ.keys():
    SHIMURADIR = os.environ['SHIMURADIR'];
    gp("SHIMURADIR=\""+SHIMURADIR+"\";");
    # print("SHIMURADIR = %s" % SHIMURADIR);
else:
    print("NO SHIMURADIR.");