
#!/bin/env bash
#

PERLBREW='http://install.perlbrew.pl'
BLASTLINUX='ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/ncbi-blast-2.6.0+-x64-linux.tar.gz'
BLASTMAC='ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/ncbi-blast-2.6.0+-x64-macosx.tar.gz'
BLASTFOLDER=blast

# PerlBrew needs to be installed to manage isolated perl environemnts if missing
command -v perlbrew >/dev/null 2>&1 || {
    echo 'Installing Perl brew...'
    curl -L ${PERLBREW} | bash
    source ~/perl5/perlbrew/etc/bashrc
    echo 'source ~/perl5/perlbrew/etc/bashrc' >> ~/.bash_profile

    perlbrew init

    echo "Do you want to install a local perl? [Y]/[N]"
    read answer
    if  [ $answer == 'Y' ]; then
        echo 'Installing perl-5.23.0';
        perlbrew install perl-5.23.0
        perlbrew use perl-5.23.0
    else
        echo 'Local perl will be used';
    fi
}

echo "Do you want to install a cpanmin locally through perlbrew? [Y]/[N]"
read answer
if  [ $answer == 'Y' ]; then
    echo 'Installing perlbrew install-cpanm';
    perlbrew install-cpanm
fi

echo "Do you want to install a cpanmin as sudo[Y]/[N]"
read answer
if  [ $answer == 'Y' ]; then
    echo 'Installing cpanmin as sudo';
  curl -L https://cpanmin.us | perl - --sudo App::cpanminus
fi

# Installing NCBI Blast tools if missing
command -v blastall >/dev/null 2>&1 || {
    echo 'Installing Blast tools...'
    curl ${BLASTMAC} -o ${BLASTFOLDER}.tar.gz
    tar -zxvf ${BLASTFOLDER}.tar.gz
}
