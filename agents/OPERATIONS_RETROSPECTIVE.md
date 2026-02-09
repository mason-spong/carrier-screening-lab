# Operations Retrospective: Execution Context Pitfalls

## What went wrong

1. We treated agent-launched background processes as if they were guaranteed to live in the same runtime context as the user's desktop/browser session.
2. We assumed `127.0.0.1` seen by the agent would always be reachable from the user's browser.
3. We relied on exact process-match strings (`pgrep`) that were too strict (relative vs absolute path differences).
4. We accepted partial BWA index state (`.bwt` present) as if indexing were complete; missing `.sa` repeatedly blocked downstream phases.
5. We used PID files that could become stale without validating liveness before trusting status.

## First-principles analysis

- Localhost is namespace-local, not globally shared.
  - `127.0.0.1` refers to "this process namespace". If agent and user app are not in the same network/runtime context, a service can be reachable from one and refused from the other.
- Process lifetime is tied to owning session semantics.
  - A background process started from one execution context can be reaped when that context ends, even if `nohup` was used.
- File presence does not imply workflow completeness.
  - For multi-artifact indexes, one file existing is insufficient proof. Completion must be defined as all required artifacts present and non-empty.
- Monitoring is only reliable if detection is robust to command-line variants.
  - Process detection must account for relative/absolute paths and wrapper commands.

## Corrective actions implemented

1. Full BWA index completeness check added before full run:
   - Require `.amb`, `.ann`, `.bwt`, `.pac`, and `.sa`.
2. UI utility split into a user-local launcher:
   - `scripts/run_status_ui_maximalist_local.sh` starts UI in the user's own terminal/session.
3. Status server hardened:
   - `scripts/serve_status_ui_maximalist.py` supports broader process-pattern matching and dynamic rendering.
4. Orchestration pattern clarified:
   - Run supervisor from user terminal with `nohup caffeinate ...` as process of record.
5. Logs elevated as source of truth:
   - `logs/mason_orchestrator.log` and `logs/mason_local_runner.log`.

## Future-safe operating rules

1. Any user-facing web UI must be launched from the user's own terminal session.
2. Any long-running analysis process of record must be started from the user's terminal, not agent-only background context.
3. Never infer completion from single index artifacts; require full artifact set.
4. Every PID file check must be paired with `ps` validation.
5. Status pages should tolerate both relative and absolute process command patterns.

## Standard launch sequence

1. Start analysis:
   - `nohup caffeinate -dimsu bash scripts/run_mason_orchestrator.sh > logs/mason_local_runner.log 2>&1 & echo $! > logs/mason_local_runner.pid`
2. Start UI:
   - `./scripts/run_status_ui_maximalist_local.sh 8788`
3. Verify:
   - `ps -p $(cat logs/mason_local_runner.pid) -o pid,etime,command` (if PID file used)
   - `tail -n 30 logs/mason_orchestrator.log`
   - Open `http://127.0.0.1:8788/`
