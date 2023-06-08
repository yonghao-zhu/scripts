#! /usr/bin/python3
# -*- conding=UTF-8 -*-
#  .--,       .--,
# ( (  \.---./  ) )
#  '.__/o   o\__.'
#     {=  ^  =}
#      >  -  <
#     /  Zhu  \
#    //  Yong \\
#   //|  Hao  |\\
#   "'\       /'"_.-~^`'-.
#      \  _  /--'         `
#    ___)( )(___
#   (((__) (__))) 

'''
    calculate three-order ifc
    run me at a Linux system
'''

import os

#------------------------------------------------------------------------
def allocate(prefix, num_jobs):

    if not os.path.isdir('./force_run'):
        os.system("mkdir ./force_run/")
    
    if num_jobs > 999:
        exit("Ary you sure? Jobs are MORE than 999? Exiting...")
    
    if num_jobs < 9:
        order = 1; order_ = "%d"
    elif (num_jobs > 9) and (num_jobs < 99):
        order = 2; order_ = "%02d"
    else:
        order = 3; order_ = "%03d"
    print("The order is %s!" %order)

    for ijob in range(1, num_jobs+1):
        jobname = "./force_run/" + prefix + "%d"%ijob
        scfname = "./force_run/" + prefix + "%d"%ijob + "/scf.in"
        
        if not os.path.isdir(jobname):
            os.system("mkdir %s" %jobname)
        
        if os.path.isfile(prefix + order_ %ijob):
            print(jobname)
            os.system("mv %s %s" %(prefix + order_ %ijob, scfname))

#------------------------------------------------------------------------
def copyout(prefix, num_jobs):

    if not os.path.isdir('./force_run'):
        exit("No ./force_run/! Exiting...")

    if not os.path.isdir('./out'):
        os.system("mkdir ./out")

    if num_jobs < 9:
        order = 1; order_ = "%d"
    elif (num_jobs > 9) and (num_jobs < 99):
        order = 2; order_ = "%02d"
    else:
        order = 3; order_ = "%03d"
    print("The order is %s!" %order)

    for ijob in range(1, 1+num_jobs):
        jobname  = "./force_run/" + prefix + "%d"%ijob
        filename = "./out/out_" + order_ %ijob + ".log"

        if not os.path.isfile(jobname + "/out.log"):
            exit("No %s/out.log! Exiting..." %jobname)
        
        print("cp %s/out.log %s" %(jobname, filename))
        os.system("cp %s/out.log %s" %(jobname, filename))

#------------------------------------------------------------------------
def checkconvergence(prefix, num_jobs):

    for i in range(num_jobs):
        ifile_name = './out/out_%03d.log' %(i+1)

        if not os.path.isfile(ifile_name):
            exit('No %s! Exitting...' %ifile_name)

        with open(ifile_name, 'r') as f:
            ifile = f.readlines()
            for iline in ifile:
                if 'convergence NOT achieved' in iline:
                    print('Not Convergence-->%s' %(i+1))


      
#------------------------------------------------------------------------
def main():
    
    prefix   = "DISP.mos2_sc.in." # DISP.BAs_sc.in.xxx 
    num_jobs = 176 # starts 1

    # 1st
    #allocate(prefix, num_jobs)

    # 2nd
    #copyout(prefix, num_jobs)

    # 3rd
    checkconvergence(prefix, num_jobs)

if __name__ == "__main__":
    main()
