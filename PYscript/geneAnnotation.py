import sys
import gffutils
db = gffutils.create_db('/Users/homere/Desktop/IGV/ostreaedulisgfffile.gff', dbfn='/Users/homere/Desktop/IGV/ostreaedulisgfffile.db', force=True)

gene = db.region(region='your_chromosome:your_position-your_position', featuretype='gene')

print(gene.attributes['gene'][0])
