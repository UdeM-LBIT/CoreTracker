mtgenes = {

    'nad1': set(['nad1', 'nd1', 'mtnd1', 'mt-nd1', 'nadh1']),
    'nad2': set(['nad2', 'nd2', 'mtnd2', 'mt-nd2', 'nadh2']),
    'nad3': set(['nad3', 'nd3', 'mtnd3', 'mt-nd3', 'nadh3']),
    'nad4': set(['nad4', 'nd4', 'mtnd4', 'mt-nd4', 'nadh4', 'lhon']),
    'nad4L': set(['nad4l', 'nd4l', 'mtnd4l', 'mt-nd4l', 'nadh4l']),
    'nad5': set(['nad5', 'nd5', 'mtnd5', 'mt-nd5', 'nadh5']),
    'nad6': set(['nad6', 'nd6', 'mtnd6', 'mt-nd6', 'nadh6']),
    'nad7': set(['nad7', 'nd7', 'mtnd7', 'mt-nd7', 'nadh7']),
    'nad8': set(['nad8', 'nd8', 'mtnd8', 'mt-nd8', 'nadh8']),
    'nad9': set(['nad9', 'nd9', 'mtnd9', 'mt-nd9', 'nadh9']),
    'nad10': set(['nad10', 'nd10', 'mtnd10', 'mt-nd10', 'nadh10']),
    'nad11': set(['nad11', 'nd11', 'mtnd11', 'mt-nd11', 'nadh11']),

    'atp1': set(['atp1', 'mtatp1', 'mt-atp1', 'atpase-1']),
    'atp3': set(['atp3', 'mtatp3', 'mt-atp3', 'atpase-3']),
    'atp4': set(['atp4', 'mtatp4', 'mt-atp4', 'atpase-4']),
    'atp6': set(['atp6', 'mtatp6', 'mt-atp6', 'atpase-6', 'su6m', 'rp']),
    'atp8': set(['atp8', 'mtatp8', 'mt-atp8', 'atpase-8', 'a6l']),
    'atp9': set(['atp9', 'mtatp9', 'mt-atp9', 'atpase-9']),

    'cox1': set(['cox1', 'coi', 'cox-1', 'mt-co1', 'mtco1', 'mtcox1', 'mt-cox1']),
    'cox2': set(['cox2', 'coii', 'cox-2', 'mt-co2', 'mtco2', 'mtcox2', 'mt-cox2']),
    'cox3': set(['cox3', 'coiii', 'cox-3', 'mt-co3', 'mtco3', 'mtcox3', 'mt-cox3']),
    'cox11': set(['cox11', 'cox-11', 'mtcox11', 'mt-cox11']),

    'ccma': set(['ccma', 'yejw']),
    'ccmb': set(['ccmb', 'yejv']),
    'ccmc': set(['ccmc', 'yeju']),
    'ccmf': set(['ccmf', 'yejr']),

    'tatc': set(['tatc', 'mttb']),
    'tata': set(['tata', 'mtta']),

    'cob': set(['mt-cyb', 'mtcyb', 'cob', 'cytb', 'uqcr3'])

}

revmtgenes = dict((v, k) for k in mtgenes.keys() for v in mtgenes[k])
try:
    import os
    curdir = os.path.dirname(os.path.realpath(__file__))
    tmpdir = os.path.join(curdir, 'templates')
    with open(os.path.join(tmpdir, "mtcorr.txt")) as corr:
        for line in corr:
            key, value = line.strip().split(None, 1)
            correspondances = value.split('|')
            key = key.lower()
            if key in revmtgenes.keys():
                mtgenes[revmtgenes[key]].update(correspondances)
            else:
                reverse_find = False
                for c in correspondances:
                    if c in mtgenes.keys():
                        reverse_find = True
                        mtgenes[c].add(key)
                if not reverse_find:
                    mtgenes[key] = set(correspondances)
except:
    print "Error parsing mtgenes alias"
    pass
revmtgenes = dict((v, k) for k in mtgenes.keys() for v in mtgenes[k])
