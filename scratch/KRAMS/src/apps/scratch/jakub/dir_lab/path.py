import os

home_dir = os.environ['HOME']
print os.path.basename(home_dir)
base_name = os.path.basename(os.getcwd())
print home_dir

data_dir = home_dir+'/simdata'
if not os.path.exists(data_dir):
    os.mkdir(data_dir)
    print "simdata directory created"
    
if not os.path.exists(data_dir +'/'+ base_name):
    os.mkdir(data_dir +'/'+ base_name)
    print base_name," directory created"