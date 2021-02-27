import math

class body:
    def __init__(self, name, mu, sma, e, IP, i, r):
        self.name = name
        self.mu = mu
        self.sma = sma
        self.e = e
        self.IP = IP
        self.i = i
        self.r = r

    def get_info(self):
        return self

### SUN ###
sun = body(
name = 'Sun',
mu = 132712440017.99,
sma = float('nan'),
e = float('nan'),
IP = float('nan'),
i = float('nan'),
r = 695990
)

### MOON ###
moon = body(
name = Moon,
mu = 4902.8005821478,
sma = 384400,
e = 0.0554,
IP = 2360592,
i = math.radians(5.16),
r = 1738.2
)

### MERCURY ###
mercury = body(
name = Mercury,
mu = 22032.080486418,
sma =57909101,
e = 0.20563661,
IP = 7600537,
i = math.radians(7.00497902),
r = 2439.7
)

### VENUS ###
venus = body(
name = Venus,
mu = 324858.59882646,
sma = 108207284,
e = 0.00676399,
IP = 19413722,
i = math.radians(3.39465605),
r = 6051.9
)

### EARTH ###
earth = body(
name = Earth,
mu = 398600.4415,
sma = 149597898,
e = 0.01673163,
IP = 31558205,
i = math.radians(0.00001531),
r = 6378.1363
)

### MARS ###
mars = body(
name = Mars,
mu = 42828.314258067,
sma = 227944135,
e = 0.09336511,
IP = 59356281,
i = math.radians(1.84969142),
r = 3397
)

### JUPITER ###
jupiter = body(
name = Jupiter,
mu = 126712767.8578,
sma = 778279959,
e = 0.04853590,
IP = 374479305,
i = math.radians(1.30439695),
r = 71492
)

### SATURN ###
saturn = body(
name = Saturn,
mu = 37940626.061137,
sma = 1427387908,
e = 0.05550825,
IP = 930115906,
i = math.radians(2.48599187),
r = 60268
)

### URANUS ###
uranus = body(
name = Uranus,
mu = 5794549.007019,
sma = 2870480873,
e = 0.04685740,
IP = 2652503938,
i = math.radians(0.77263783),
r = 25559
)

### NEPTUNE ###
neptune = body(
name = Neptune,
mu = 6836534.0638793,
sma = 4498337290,
e = 0.00895439,
IP = 5203578080,
i = math.radians(1.77004347),
r = 25269
)

### PLUTO ###
pluto = body(
name = Pluto,
mu = 981.600887707,
sma = 5907150229,
e = 0.24885238,
IP = 7830528509,
i = math.radians(17.14001206),
r = 1162
)
