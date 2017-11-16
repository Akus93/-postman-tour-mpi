# Skrypt konwertuje pliki z https://snap.stanford.edu/data/#road na pasujace do programu
# Plik nazwac input.txt, wyjscie to output.txt

verticies = {}
result = []

with open('input.txt', 'r') as file:
    for line in file:
        row = tuple(map(int, line.split()))
        try:
            verticies[row[0]].append(row[1])
        except KeyError:
            verticies[row[0]] = []
            verticies[row[0]].append(row[1])

for key, value in verticies.items():
    for i in value:
        verticies[i].remove(key)

for key, value in verticies.items():
    for i in value:
        result.append((key, i))

result.sort(key=lambda item: item[0])

with open('output.txt', 'w') as output:
    output.writelines(['{} {}\n'.format(*x) for x in result])
