import sys,re
# usage: ctx.py file pattern [width]
f,p=sys.argv[1],sys.argv[2]; w=int(sys.argv[3]) if len(sys.argv)>3 else 250
t=open(f).read()
for m in re.finditer(p,t):
    ln=t.count('\n',0,m.start())+1
    print(f'--- L{ln}:', t[max(0,m.start()-w):m.end()+w].replace('\n',' | '))
