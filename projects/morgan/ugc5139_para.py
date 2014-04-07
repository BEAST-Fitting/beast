import sys

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('must have a chunk number')
        sys.exit(1)
    chunk = int(sys.argv[1])
    import ugc5139
    g = '{project:s}/{project:s}_seds.grid.hd5'.format(project=ugc5139.project)
    ugc5139.run_chunk_fit(ugc5139.project, g, chunk)
