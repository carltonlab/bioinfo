import json
import argparse

# you might want to preprocess the clustal file with the following to get each 
# record on one line:
# perl -ane '$h{$F[0]}.=$F[1];if(eof){for $key(keys %h) {print $key," "x(40-(length($key))),$h{$key},"\n"}};' ORIGINAL.clu | grep -v CLUSTAL > ONELINEEACH.clu

def getjson(jsonfilename):
    jsonfile=open(jsonfilename)
    dat=json.load(jsonfile);jsonfile.close()
    q8={x['id'][5:35]:x['q8'] for x in dat}
    disorder={x['id'][5:35]:x['disorder'] for x in dat}
    return(q8,disorder)

def cluprint(clufilename,q8,disorder):
    darkness=50 #how dark the "disorder" color background of text should be relative to 255=total
    print('''
<!DOCTYPE html>
<html>
<head>

<style>
.G { color: red !important; }
.H { color: orange !important; }
.I { color: brown !important; }
.T { color: green !important; }
.E { color: blue !important; }
.B { color: purple !important; }
.S { color: darkturquoise !important; }
.C { color: fuchsia !important; }
</style>

</head>

<body>

<h1>Clustal-Netsurf</h1>
<pre>
Color code (for q8):
<span class=G>Red: G (3-turn helix)</span>
<span class=H>Orange: H (4-turn helix)</span>
<span class=I>Brown: I (5-turn helix)</span>
<span class=T>Green: T (hydrogen-bonded turn)</span>
<span class=E>Blue: E (extended strand in parallel and/or antiparallel beta-sheet)</span>
<span class=B>Purple: B (residue in isolated beta-bridge)</span>
<span class=S>Darkturquoise: S (bend)</span>
<span class=C>Magenta: C (coil)</span>

- Magenta background strength correlates with higher disorder parameter

    ''')
    count={x:0 for x in q8.keys()}
    clufile=open(clufilename)
    clu=list(clufile);clufile.close()
    for line in clu:
        splitline=line.split()
        try:
            name=splitline[0].replace(".","_")
        except:
            name=""
        if name in q8.keys():
            spc=" "*(40-len(name))
            print(name,spc,end="")
            for char in splitline[1]:
                if(char != "-"):
                    bval="{:02x}".format(255-int(darkness*disorder[name][count[name]]))
                    bval="#FF"+(bval+"FF") #subtract from green depending on disorder
                    print('<span style="background-color: ',bval,'" class="',q8[name][count[name]],'">',char,'</span>',sep="",end="")
                    count[name]=count[name]+1
                else:
                    print(char,end="")
            print()
#        else:
#            print("xx_xx",line)
            #print("    ",line)
    print('''
    </pre>
    </body>
    </html>
    ''')


def main():
    parser=argparse.ArgumentParser(description='Output a color-coded MSA file given clustal and json')
    parser.add_argument('clufilename',type=str, help='textfile containing a Clustal-format MSA')
    parser.add_argument('jsonfilename',type=str, help='textfile containing a JSON output of netsurfp2 of proteins in MSA')
    args=parser.parse_args()
    q8,disorder=getjson(args.jsonfilename)
    cluprint(args.clufilename,q8,disorder)


if __name__=="__main__":
    main()
