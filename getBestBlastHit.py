#!/usr/bin/python
import argparse

def getBestHits(file, out_file):
    ## creating output file
    outfile = open(out_file, 'w')
    query = open(file, 'r')
    outfile.write('query' + '\t' + 'hit' + '\t' + 'evalue' + '\n')

    id_list = []
    uniqueID = ''
    for line in query:
        line = line.strip().split()
        queryID = line[0]
        hit = line[1]
        evalue = line[9]
        #bit_score = line[len(line)-1]
        #print hit
        if  queryID != uniqueID:
            #print queryID + '\t' + uniqueID
            outfile.write(queryID + '\t' + hit + '\t' + evalue + '\n')
            uniqueID = queryID
            #print uniqueID
            id_list.append(uniqueID)
    query.close()
    outfile.close()
    return id_list



def main():
    parser = argparse.ArgumentParser(description='Get Best Blast Hit \
    using ncbi-blast outpu in tabular format (6)')
    parser.add_argument('--query', required=True, metavar='tab', nargs='+',
                        help='[REQUIRED] Blast output file in tab format (6)')
    parser.add_argument('-o', metavar='output', type=str,
                        help='[optional] extension of output file(s). if not specified,\
                         the output file(s) name(s) will be the same as query file(s)\
                          with "_BestHit.txt" extension ')
    args = parser.parse_args()

    for query in args.query:
        if args.o:
            outFileName = query + args.o
        else:
            outFileName = query + '_BestHit.txt'

        query_IDs = getBestHits(query, outFileName)

        print str(len(query_IDs)) + ' genes with blast hits'

if __name__ == "__main__":
    main()
