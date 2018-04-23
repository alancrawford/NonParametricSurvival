struct np_survival
	t :: Vector
	d :: Vector
	m :: Vector
	r :: Vector
end 

struct survivor
	surv  :: Vector
	V  	  :: Vector
	se 	  :: Vector
	sigma :: Vector
end 

struct cumu_hazard
	ch :: Vector
	V  :: Vector
	se :: Vector
	sigma :: Vector
end 