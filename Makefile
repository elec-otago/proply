# Build every propeller in props/ and render it to images/<name>.png.
#
#   make gallery            design all props and render each to images/
#   make steps              design all props (STEP files in build/out/)
#   make summaries          design all props (YAML summaries in build/out/)
#   make wasm               build the WebAssembly package (proply-rs/pkg/)
#   make build-date         stamp the web build label (web/build.js)
#   make check-wasm         type-check the lib for wasm32 without wasm-pack
#   make clean              remove generated STEP, YAML and PNG files
#
# Designs use the coupled lifting-line solver by default.  Switch design
# modes by overriding DESIGN_FLAGS, e.g.
#   make gallery DESIGN_FLAGS="--naca --bem --n 40 --element-count 30"

DESIGN_FLAGS ?= --naca --lifting-line --mech-thickness --n 40 --element-count 30

# Never delete intermediate artefacts (the .ply meshes in particular are
# the render inputs and each design run is expensive).
.SECONDARY:

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
# match (and ntm_28_26_1200Kv.json omits it).  --mesh-file writes the
# triangle mesh (PLY) the headless render-step previewer draws.
build/out/%.step build/out/%.yml build/out/%.ply &: props/%.json $(STAMP)
	@mkdir -p $(dir $@)
	cargo run --release -p proply-rs -- $(DESIGN_FLAGS) --log build/out/$*.log --step-file=build/out/$*.step --mesh-file=build/out/$*.ply --param=$<

# The headless render-step utility rasterises the design's triangle mesh
# (written beside the STEP) to the gallery PNG — no CAD kernel, no display.
images/%.png: build/out/%.step build/out/%.ply build/out/%.yml
	@mkdir -p images
	cargo run --quiet --release -p render-step -- --step $< --png $@
	test -f $@

clean:
	rm -f $(STEPS) $(YAMLS) $(PNGS) $(STAMP)

# WebAssembly: build the browser package into proply-rs/web/pkg/ (so the
# web/ directory is a self-contained static site), or just type-check the
# lib for the wasm32 target.
WASM_TARGET := wasm32-unknown-unknown

wasm:
	wasm-pack build proply-rs --target web --release
	rm -rf proply-rs/web/pkg
	mv proply-rs/pkg proply-rs/web/pkg
	rm -f proply-rs/web/pkg/.gitignore  # wasm-pack ignores its own output;
	                                    # the web/pkg copy is meant to be committed

# Build stamp for the web demo: the deployed page shows "build
# yyyy-mm-dd.xx" — the last commit's date (`%cs`) and its per-day build
# number (the count of commits on that date, so each build on a day
# increments xx and a new day starts at .01).  Written to
# proply-rs/web/build.js, which main.js renders; regenerate before
# committing web changes so the label matches the deployed sources.
build-date:
	@DATE=$$(git log -1 --format=%cs); \
	COUNT=$$(git log --format=%cs | grep -cx "$$DATE"); \
	XX=$$(printf '%02d' "$$COUNT"); \
	{ printf '// Build stamp: "yyyy-mm-dd.xx" — the last commit date and its per-day\n// build number (commits on that date).  Regenerate with make build-date\n// before committing web changes.\n'; \
	  printf 'export const BUILD = "%s.%s";\n' "$$DATE" "$$XX"; \
	} > proply-rs/web/build.js; \
	echo "web build stamp: $$DATE.$$XX -> proply-rs/web/build.js"

check-wasm:
	cargo check -p proply-rs --target $(WASM_TARGET)

.PHONY: all steps summaries gallery clean wasm check-wasm build-date
.DELETE_ON_ERROR:
