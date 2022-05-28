def load_conditions(filename):
    file = open(filename, 'r')
    lines = file.readlines()
    arr = []
    print(filename)
    print(len(lines))
    for one_line in lines:
        line = one_line.strip()
        if line == '0':
            arr.append([1,0,0])
        if line == '1':
            arr.append([0,1,0])
        if line == '2':
            arr.append([0,0,1])
    print(len(arr))
    return arr
    