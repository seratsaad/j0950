# Download data for proprietary HST data
# copied script from this site: https://archive.stsci.edu/ftp.html
# Updated to put the data in a specific directory
# First go to e.g., a website like this: https://archive.stsci.edu/proposal_search.php?mission=hst&id=16765
# And request the data be sent to the staging area
#
# J. Runnoe
# 04/13/2022
# Written for project 16765 (J0950 UV spec)

# import block
import ftplib


if __name__=="__main__":

    # set specific user and data information
    # password is for MySTScI website account
    user     = 'jessie.c.runnoe@vanderbilt.edu'
    password = ''
    stagedir = '/stage/runnoejc/runnoejc86339' 
    datadir  = '../fits/'

    ftps = ftplib.FTP_TLS('archive.stsci.edu')
    ftps.login(user=user, passwd=password)
    ftps.prot_p() # This is a really good idea :)
    ftps.cwd('stage')
    ftps.cwd(stagedir) # stagedir is something like 'anonymous/anonyumous12345'

    filenames = ftps.nlst()
    for filename in filenames:
        print("getting " + filename)
        with open(datadir+filename, 'wb') as fp: 
            ftps.retrbinary('RETR {}'.format(filename), fp.write)

