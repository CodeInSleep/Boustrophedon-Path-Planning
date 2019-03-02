epsilon = 1e-4


BEGINEVENT = 'BEGIN EVENT'
ENDEVENT = 'END EVENT'
INEVENT = 'IN EVENT'
OUTEVENT = 'OUT EVENT'
CEILINGEVENT = 'CEILING EVENT'
FLOOREVENT = 'FLOOR EVENT'
INOUTEVENTS = [INEVENT, OUTEVENT]
MIDEVENTS = [CEILINGEVENT, FLOOREVENT]

def approxEqual(val1, val2):
	return abs(val1 - val2) <= ((abs(val1) if abs(val2) < abs(val1)
		else abs(val2))*epsilon)