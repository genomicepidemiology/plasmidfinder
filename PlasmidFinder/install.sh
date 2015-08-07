# Installing cpanm if missing
command -v cpanm >/dev/null 2>&1 || {
    if [[ "${OSTYPE}" == 'linux'* ]]; then
            wget -O - https://cpanmin.us | perl - --sudo App::cpanminus
    else
            curl -L https://cpanmin.us | perl - --sudo App::cpanminus
    fi
}

# Installing BioPerl if missing
perl -MBio::Seq -e 0 >/dev/null 2>&1 || {
    cpanm -nf BioPerl
    #cpanm http://search.cpan.org/CPAN/authors/id/C/CJ/CJFIELDS/BioPerl-1.6.901.tar.gz
}

perl Makefile.PL
make
make install
make clean

# Installing NCBI Blast tools if missing
BLASTLINUX='ftp://ftp.ncbi.nlm.nih.gov/blast/executables/release/LATEST/blast-2.2.26-x64-linux.tar.gz'
BLASTMAC='ftp://ftp.ncbi.nlm.nih.gov/blast/executables/release/LATEST/blast-2.2.26-universal-macosx.tar.gz'
BLASTFOLDER=blast
command -v blastall >/dev/null 2>&1 || {
    echo 'Installing Blast tools...'

    if [[ "${OSTYPE}" == 'linux'* ]]; then
        wget ${BLASTLINUX}
    else
        # TODO Include versions for all OS: BSD, etc...
        curl ${BLASTMAC} > ${BLASTFOLDER}.tar.gz
        tar -zxvf ${BLASTFOLDER}.tar.gz
    fi

}