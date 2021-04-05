#3d Lennard-Jones melt

variable        lj_density equal 0.7
variable        lj_temp equal 11.0

log             melt.d${lj_density}.T${lj_temp}.txt
# Единицы измерения (по умолчанию LJ)
units		lj

# Сколько атрибутов есть у атомов
atom_style	atomic # нет зарядов, нет молекул, у атомов есть только тип и масса
# другие atom_style:
# * charge (есть заряды)
# * molecular (есть молекулы)
# ...
# * full (есть все параметры, которые LAMMPS поддерживает)

# тип и характерный размер кристаллической решётки
lattice		fcc ${lj_density}
# * fcc - face-centered cubic
# * bcc - bulk-centered cubic
# * sc - simple cubic
# etc

# многие команды считают параметр кристаллической решётки
# характерным масштабом длины

# units lj: параметр решётки - это плотность
# units !lj: параметр решётки - это длина в единицах длины

# region - обозначает область пространства идентификатором
region		simbox block 0 10 0 10 0 10

# box - идентификатор
# block - форма пр. параллелепипеда (может быть sphere, cylinder, ...)
# 0 10 0 10 0 10 - границы блока (ниж. верх. ниж. верх. ниж. верх.)
# не указан параметр по умолчанию - units lattice
# units lattice - числа в единицах размера решётки
# units box - числа в единицах длины

# создание расчётной ячейки внутри ранее отмеченного региона
create_box	1 simbox

# 1 - число типов атомов в расчёте
# simbox - имя региона

# заполнить расчётную область атомами
create_atoms	1 box

# 1 - тип создаваемых атомов
# box - создать атомы в узлах решётки, находящихся внутри расчётной области
# на месте box может быть
# * region <имя> - создать атомы в узлах решётки, находящихся внутри региона <имя>
# * random <N> <seed> <region-ID>
# и т.д. - см. документацию

# задать массы для типов атомов
mass		1 1.0 # для типа 1 установлена масса 1.0

# если не задать массы - при запуске расчёта будет ошибка all masses are not set


velocity	all create 3.0 87287 # задать скорости группе all случайно, средняя кинетическая энергия соответствует температуре 3.0

# тип потенциала
pair_style	lj/cut 2.5 # потенциал LJ с обрезанием, радиус обрезания 2.5
pair_coeff	1 1 1.0 1.0 2.5

# pair_coeff <type1> <type2> <pair_parameters...>
# в <type1> или <type2> могут быть *

# параметры составления списка соседей
neighbor	0.3 bin # радиус списка соседей на 0.3 больше радиуса обрезания
neigh_modify	every 20 delay 0 check no # как часто перестраивать список соседей

# fix - изменять координаты, скорости и силы некоторым образом (обычно)
fix		1 all nve

# fix <fix-id> <group-id> <fix-style> <fix-parameters...>
# fix-id: 1
# group-id: all
# fix-style: nve
# fix nve не требует параметров

#fix thermostat all temp/berendsen ${lj_temp} ${lj_temp} 0.2

# compute - в ходе расчёта вычислять некоторую функцию от свойств атомов



# выводить термодинамические параметры на шагах с номером, кратным 50
#compute 1 all pe
#thermo	1000
dump dump_1 all custom 10 yt.dump id type xu yu zu vx vy vz
#dump dump_2 all custom 10 a.dump id type pe


# пропущено: thermo_style (по умолчанию - )
timestep 0.005
# пропущено: run_style (по умолчанию - verlet)

# провести расчёт длиной 10000 шагов
run		10000

# сбросить номер шага по времени
reset_timestep  0
run             10000


#dump		id all atom 50 dump.melt

#dump		2 all image 25 image.*.jpg type type &
#		axes yes 0.8 0.02 view 60 -30
#dump_modify	2 pad 3

#dump		3 all movie 25 movie.mpg type type &
#		axes yes 0.8 0.02 view 60 -30
#dump_modify	3 pad 3

