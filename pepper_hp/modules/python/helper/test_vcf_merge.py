def test_vcf(filename):
    with open(filename, encoding='utf-8') as handle:
        for index, line in enumerate(handle):
            line = line.replace('\n', '')
            # Already read meta and header in self.__init__
            if line.startswith('#'):
                continue
            list_line = line.split("\t")
            ref = list_line[3]
            alt = list_line[4]
            alts = alt.split(',')
            min_len = min([len(alt) for alt in alts])
            print(len(ref), min_len)
            if len(ref) == 1 or min_len == 1:
                continue
            else:
                print(ref, alt)
                print("ERROR", line)

            # if len(ref) > 50 or len(alt) > 50:
            #     print(ref, alt)
            #     print("SV", line)