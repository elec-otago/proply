/* tslint:disable */
/* eslint-disable */

/**
 * One finished design: the STEP (AP242) document, the YAML summary and
 * the headline numbers.
 */
export class DesignOutput {
    private constructor();
    free(): void;
    [Symbol.dispose](): void;
    power: number;
    rpm: number;
    step: string;
    thrust: number;
    torque: number;
    /**
     * Empty when the design reached its operating point; otherwise an
     * explicit note describing the closest achievable design.
     */
    warning: string;
    yaml: string;
}

/**
 * A design session with a warm polar cache, kept across design calls.
 */
export class PropSession {
    free(): void;
    [Symbol.dispose](): void;
    /**
     * The whole cache as a JSON document (export/migration escape hatch).
     */
    cache_to_json(): string;
    /**
     * Run one full design from JSON design parameters (the same format
     * and validation as the CLI's `--param` file).
     */
    design(params_json: string): DesignOutput;
    /**
     * Insert one pre-existing polar record — the bulk startup hydration
     * path.  Records inserted here are not reported by
     * [`PropSession::take_new_json`] (the host already has them).
     */
    hydrate_entry(key: string, alpha: Float64Array, cl: Float64Array, cd: Float64Array): void;
    /**
     * Hydrate from a full cache document ([`PropSession::cache_to_json`]
     * format), replacing any current contents.
     */
    hydrate_json(json: string): void;
    constructor();
    /**
     * Number of cached polars (e.g. after startup hydration).
     */
    polar_count(): number;
    /**
     * Install the host's per-polar persistence hook: called synchronously
     * for every polar the moment it is freshly calculated — a good sweep
     * or a degenerate failure marker — with the cache key and the
     * (alpha, cl, cd) arrays.  The host writes each record to its
     * IndexedDB cache immediately, so a design interrupted mid-way keeps
     * every completed sweep, exactly like the native CLI's per-polar disk
     * checkpoint.  The hook replaces any previously installed one.
     */
    set_on_polar(on_polar: Function): void;
    /**
     * A JSON map `key -> {alpha, cl, cd}` of the polars simulated since
     * the last call — what the host should persist (e.g. into IndexedDB).
     * The session keeps every polar for future designs.
     */
    take_new_json(): string;
}

export type InitInput = RequestInfo | URL | Response | BufferSource | WebAssembly.Module;

export interface InitOutput {
    readonly memory: WebAssembly.Memory;
    readonly __wbg_designoutput_free: (a: number, b: number) => void;
    readonly __wbg_get_designoutput_power: (a: number) => number;
    readonly __wbg_get_designoutput_rpm: (a: number) => number;
    readonly __wbg_get_designoutput_step: (a: number) => [number, number];
    readonly __wbg_get_designoutput_thrust: (a: number) => number;
    readonly __wbg_get_designoutput_torque: (a: number) => number;
    readonly __wbg_get_designoutput_warning: (a: number) => [number, number];
    readonly __wbg_get_designoutput_yaml: (a: number) => [number, number];
    readonly __wbg_propsession_free: (a: number, b: number) => void;
    readonly __wbg_set_designoutput_power: (a: number, b: number) => void;
    readonly __wbg_set_designoutput_rpm: (a: number, b: number) => void;
    readonly __wbg_set_designoutput_step: (a: number, b: number, c: number) => void;
    readonly __wbg_set_designoutput_thrust: (a: number, b: number) => void;
    readonly __wbg_set_designoutput_torque: (a: number, b: number) => void;
    readonly __wbg_set_designoutput_warning: (a: number, b: number, c: number) => void;
    readonly __wbg_set_designoutput_yaml: (a: number, b: number, c: number) => void;
    readonly propsession_cache_to_json: (a: number) => [number, number];
    readonly propsession_design: (a: number, b: number, c: number) => [number, number, number];
    readonly propsession_hydrate_entry: (a: number, b: number, c: number, d: any, e: any, f: any) => void;
    readonly propsession_hydrate_json: (a: number, b: number, c: number) => void;
    readonly propsession_new: () => number;
    readonly propsession_polar_count: (a: number) => number;
    readonly propsession_set_on_polar: (a: number, b: any) => void;
    readonly propsession_take_new_json: (a: number) => [number, number];
    readonly __wbindgen_malloc: (a: number, b: number) => number;
    readonly __wbindgen_realloc: (a: number, b: number, c: number, d: number) => number;
    readonly __wbindgen_exn_store: (a: number) => void;
    readonly __externref_table_alloc: () => number;
    readonly __wbindgen_externrefs: WebAssembly.Table;
    readonly __wbindgen_free: (a: number, b: number, c: number) => void;
    readonly __externref_table_dealloc: (a: number) => void;
    readonly __wbindgen_start: () => void;
}

export type SyncInitInput = BufferSource | WebAssembly.Module;

/**
 * Instantiates the given `module`, which can either be bytes or
 * a precompiled `WebAssembly.Module`.
 *
 * @param {{ module: SyncInitInput }} module - Passing `SyncInitInput` directly is deprecated.
 *
 * @returns {InitOutput}
 */
export function initSync(module: { module: SyncInitInput } | SyncInitInput): InitOutput;

/**
 * If `module_or_path` is {RequestInfo} or {URL}, makes a request and
 * for everything else, calls `WebAssembly.instantiate` directly.
 *
 * @param {{ module_or_path: InitInput | Promise<InitInput> }} module_or_path - Passing `InitInput` directly is deprecated.
 *
 * @returns {Promise<InitOutput>}
 */
export default function __wbg_init (module_or_path?: { module_or_path: InitInput | Promise<InitInput> } | InitInput | Promise<InitInput>): Promise<InitOutput>;
