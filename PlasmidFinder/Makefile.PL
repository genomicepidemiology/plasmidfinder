use 5.010000;
use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.

WriteMakefile(
    NAME              => 'PlasmidFinder',
    VERSION_FROM      => 'lib/PlasmidFinder.pm', # finds $VERSION
    PREREQ_PM         => {
        # 'BioPerl' => 1.006901,
        'Data::Dumper' => 0,
        'File::Temp' => 0,
        'Getopt::Long' => 0,
        'Try::Tiny::Retry' => 0,
    }, # e.g., Module::Name => 1.1
    # META_MERGE => {
    #     'meta-spec' => { version => 2 },
    #      resources => {
    #          repository => {
    #              type => 'git',
    #              url  => 'https://github.com/dwimperl/Task-DWIM.git',
    #              web  => 'https://github.com/dwimperl/Task-DWIM',
    #          },
    #      },
    # },
    EXE_FILES => [ qw(lib/PlasmidFinder.pm) ],
    ($] >= 5.005 ?     ## Add these new keywords supported since 5.005
      (ABSTRACT_FROM  => 'lib/PlasmidFinder.pm', # retrieve abstract from module
       AUTHOR         => 'Jose Luis Bellod Cisneros <cisneros@cbs.dtu.dk>,
                          Rosa Allesøe <rosa@cbs.dtu.dk>') : ()),
);
