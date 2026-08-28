# Build every propeller in props/ and render it to images/<name>.png.
#
#   make gallery            design all props and render each to images/
#   make steps              design all props (STEP files in build/out/)
#   make summaries          design all props (YAML summaries in build/out/)
#   make wasm               build the WebAssembly package (proply-rs/pkg/)
#   make check-wasm         type-check the lib for wasm32 without wasm-pack
#   make clean              remove generated STEP, YAML and PNG files
#
# Designs use the coupled lifting-line solver by default.  Switch design
# modes by overriding DESIGN_FLAGS, e.g.
#   make gallery DESIGN_FLAGS="--naca --bem --n 40 --resolution 30"

DESIGN_FLAGS ?= --naca --lifting-line --n 40 --resolution 30

PROPS  := $(wildcard props/*.json)
STEPS  := $(PROPS:props/%.json=build/out/%.step)
YAMLS  := $(PROPS:props/%.json=build/out/%.yml)
PNGS   := $(PROPS:props/%.json=images/%.png)

# Changing DESIGN_FLAGS must redesign every prop: the flags are recorded in
# a stamp the design rules depend on, refreshed only when they change.
STAMP := build/out/.design_flags

all: gallery

steps: $(STEPS)

summaries: $(YAMLS)

gallery: $(PNGS)

$(STAMP): Makefile
	@mkdir -p $(dir $@)
	@printf '%s\n' "$(DESIGN_FLAGS)" > $@.tmp
	@if cmp -s $@.tmp $@; then rm -f $@.tmp; else mv $@.tmp $@; echo "design flags changed -> redesigning all props"; fi

# One design run writes both artefacts (the STEP model and the YAML
# summary), so they are a grouped target: the design is rerun whenever
# either output is missing or outdated.  --step-file pins the output name
# to the JSON file stem: the "name" field inside the JSON does not always
# match (and ntm_28_26_1200Kv.json omits it).
build/out/%.step build/out/%.yml &: props/%.json $(STAMP)
	@mkdir -p $(dir $@)
	cargo run --release -p proply-rs -- $(DESIGN_FLAGS) --step-file=build/out/$*.step --param=$<

# freecadcmd forwards script arguments only when each is preceded by --pass,
# and it crashes during Qt teardown *after* the image is saved, so the exit
# status is ignored and success is judged by the PNG existing.
images/%.png: build/out/%.step build/out/%.yml props/renderprop.py
	@mkdir -p images
	-freecadcmd props/renderprop.py --pass --step --pass $< --pass --png --pass $@
	test -f $@

clean:
	rm -f $(STEPS) $(YAMLS) $(PNGS) $(STAMP)

# WebAssembly: build the browser package into proply-rs/pkg/ (imported by
# proply-rs/web/), or just type-check the lib for the wasm32 target.
WASM_TARGET := wasm32-unknown-unknown

wasm:
	wasm-pack build proply-rs --target web --release

check-wasm:
	cargo check -p proply-rs --target $(WASM_TARGET)

.PHONY: all steps summaries gallery clean wasm check-wasm
.DELETE_ON_ERROR:
