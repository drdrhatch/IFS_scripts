def get_Ip(gfile):
    f = open(gfile,'r')
    data = f.read().split('\n')
    Ip = 0.0
    for i in data:
        if 'PLASMA' in i:
            temp = i.split()
            Ip = float(temp[-1])
            print("Ip",Ip)
    if Ip==0.0:
        print("Error, couldn't fine Ip.")
        stop
    return Ip    


