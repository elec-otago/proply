# For development purposes, use the 'develop' target to install this package from the source repository (rather than the public python packages). This means that local changes to the code will automatically be used by the package on your machine.
# To do this type
#     make develop
# in this directory.
TARGET=test_prop
RESOLUTION=30
BUILDIR=build

develop:
	sudo python3 setup.py develop

# Build using the Blade Element Momentum method
test:
	mkdir -p ${BUILDIR}
	proply --arad --bem  --n 40 --resolution ${RESOLUTION} --dir=${BUILDIR} --param='props/${TARGET}.json'
	meshlabserver -i ${BUILDIR}/${TARGET}_blade.stl -o ${BUILDIR}/${TARGET}_blade.stl -s meshclean.mlx

lint:
	pylint --extension-pkg-whitelist=numpy --ignored-modules=numpy,tart_tools,dask,dask.array --extension-pkg-whitelist=astropy --extension-pkg-whitelist=dask proply

test_upload:
	rm -rf proply.egg-info dist
	python3 setup.py sdist
	twine upload --repository testpypi dist/*

upload:
	rm -rf proply.egg-info dist
	python3 setup.py sdist
	twine upload --repository pypi dist/*

.PHONY: run build
run:
	xhost +
	docker-compose run --rm proply
build:
	docker-compose build # --no-cache

# proply --bem --naca --n 40 --resolution 40 --param /props/flywoo_robo_rb1202.5_11500kv.json --auto
