epsilon = 1e-4

def approxEqual(val1, val2):
	return abs(val1 - val2) <= ((abs(val1) if abs(val2) < abs(val1)
		else abs(val2))*epsilon)