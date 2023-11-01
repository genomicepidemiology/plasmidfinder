import ruamel.yaml

# Create an ordered dictionary for each section
package = ruamel.yaml.comments.CommentedMap()
package['name'] = 'plasmidfinder'
package['version'] = '2.1.6'

source = ruamel.yaml.comments.CommentedMap()
source['url'] = 'https://bitbucket.org/genomicepidemiology/plasmidfinder/get/2.1.6.tar.gz'

build = ruamel.yaml.comments.CommentedMap()
build['number'] = 1
build['noarch'] = 'generic'

requirements = ruamel.yaml.comments.CommentedMap()
requirements['host'] = ['python>=3.5', 'wget', 'kma']
requirements['run'] = ['python>=3.5', 'wget', 'biopython', 'tabulate', 'cgecore', 'blast']

about = ruamel.yaml.comments.CommentedMap()
about['home'] = 'https://bitbucket.org/genomicepidemiology/plasmidfinder'
about['summary'] = 'PlasmidFinder service contains one python script plasmidfinder.py which is the script of the latest version of the PlasmidFinder service. The service identifies plasmids in total or partial sequenced isolates of bacteria.'
about['license'] = 'Apache-2.0'

extra = ruamel.yaml.comments.CommentedMap()
identifiers = ruamel.yaml.comments.CommentedMap()
extra['identifiers'] = identifiers

# Create a dictionary for the entire YAML content
data = ruamel.yaml.comments.CommentedMap()
data['package'] = package
data['source'] = source
data['build'] = build
data['requirements'] = requirements
data['about'] = about
data['extra'] = extra

# Serialize the data to YAML and print it
yaml = ruamel.yaml.YAML(typ='safe', pure=True)
yaml_str = yaml.dump_all(data).replace("\"{{", "{{").replace("}}\"", "}}")
print(yaml_str)

with open('conda/meta.yaml', 'w') as f:
    f.write(yaml_str)