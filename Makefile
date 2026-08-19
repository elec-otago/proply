# Build every propeller in props/ and render it to images/<name>.png.
#
#   make gallery            design all props and render each to images/
#   make steps              design all props (STEP files in build/out/)
#   make clean              remove generated STEP and PNG files
#
# The design flags can be overridden on the command line, e.g.
#   make gallery DESIGN_FLAGS="--lifting-line --ar 6"

DESIGN_FLAGS ?= --naca --bem --n 40 --resolution 30

PROPS := $(wildcard props/*.json)
STEPS := $(PROPS:props/%.json=build/out/%.step)
PNGS  := $(PROPS:props/%.json=images/%.png)

all: gallery

steps: $(STEPS)

gallery: $(PNGS)

# --step-file pins the output name to the JSON file stem: the "name" field
# inside the JSON does not always match (and ntm_28_26_1200Kv.json omits it).
build/out/%.step: props/%.json
	@mkdir -p $(dir $@)
	cargo run --release -p proply-rs -- $(DESIGN_FLAGS) --step-file=$@ --param=$<

# freecadcmd forwards script arguments only when each is preceded by --pass,
# and it crashes during Qt teardown *after* the image is saved, so the exit
# status is ignored and success is judged by the PNG existing.
images/%.png: build/out/%.step props/renderprop.py
	@mkdir -p images
	-freecadcmd props/renderprop.py --pass --step --pass $< --pass --png --pass $@
	test -f $@

clean:
	rm -f $(STEPS) $(PNGS)

.PHONY: all steps gallery clean
.DELETE_ON_ERROR:
