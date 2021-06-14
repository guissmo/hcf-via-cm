# This file takes care of logfiles and saves history.

def dumper(oh, i):
     if i in oh.keys():
         return dumps(oh[i])
     else:
         return 0

def undumper(l, i):
     if l[i] != 0:
         return loads(l[i])
     else:
         return 0

def save_history():
    """save the IPython history to a plaintext file"""
    lines = []
    # get previous lines
    # this is only necessary because we truncate the history,
    # otherwise we chould just open with mode='a'
    if os.path.exists(histfile):
        with open(histfile, 'r') as f:
            lines = f.readlines()

    # add any new lines from this session
    lines.extend(record[2] + '\n' for record in ip.history_manager.get_range())

    with open(histfile, 'w') as f:
        # limit to LIMIT entries
        f.writelines(lines[-LIMIT:])
    
    print("Saving input history to %s" % histfile)

    with open(dumpfile, 'wb') as f:
        l = [0] + [dumper(_oh, i) for i in range(1,max(_oh)+1)];
        f.write(dumps(l));

    print("Saving output dumps to %s" % dumpfile)

def read_dump(s):
    f = open(s, 'r');
    l = loads(f.read());
    return [ undumper(l, i) for i in range(1, len(l)) ];

def save_var(s, varname=None):
    if varname == None:
        varname = str(abs(hash(s)))
    F = os.path.join(os.path.expanduser(time.strftime(LOGFILEDIR+"/"+sessionname+"-"+varname+".sagevardump")));
    with open(F, 'wb') as f:
        f.write(dumps(s));
