//! The order a flowsheet's units run in.
//!
//! **Both orders are the class's, and the default is the one without a graph.**
//! `ProcessSystem.useGraphBasedExecution` is `false` (`ProcessSystem.java:265`), and
//! `runSequential` reads it as
//!
//! ```text
//! executionOrder = useGraphBasedExecution ? getTopologicalOrder() : unitOperations
//! ```
//!
//! so a flowsheet that has not asked runs in the order its instances were **added**. The
//! topological order is reached only by a caller that sets the flag, and it is not a Kahn walk:
//! `ProcessGraph.getTopologicalOrder` is a **post-order depth-first traversal of the nodes in
//! their own order**, reversed, with `edge.isBackEdge()` edges skipped.
//!
//! **The recycle edges are removed from the ordering graph, and this refutes the plan's premise.**
//! The plan for this item said the topological walk must keep them, "a recycle is still a
//! dependency: `mix1` must run after `sep1`". The class does the opposite, and it is coherent: a
//! torn stream's value is the *previous* iteration's, so the unit reading it is not waiting for
//! this iteration's producer. That is what tearing is, and the class expresses it by marking the
//! cycle-closing edge a back edge and skipping it. A port that kept them would answer a different
//! question from the machine it ports.
//!
//! **What is not ported**: `ProcessGraph`'s Tarjan strongly-connected-component pass,
//! `partitionForParallelExecution` and the level machinery. The sequential path - the one this
//! executor is a port of - never calls them: they exist for `runOptimized` and `runDataflow`,
//! which are named in `ROADMAP.md` as not ported with those classes.

use std::collections::HashMap;

use azoth_core::{AzothError, Result};

use crate::flowsheet::Flowsheet;

/// Which of the class's two orders to run.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ExecutionOrder {
    /// The flowsheet's own declaration order. **The class's default**, because
    /// `useGraphBasedExecution` is `false`.
    #[default]
    Insertion,
    /// A depth-first topological order, reached only when a caller sets the flag.
    Topological,
}

/// The instance indices, in the order they run.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if a connection names an instance the flowsheet does not
///   declare, or if the topological walk meets a cycle that no `[[recycles]]` entry tears.
pub fn execution_order(flowsheet: &Flowsheet, mode: ExecutionOrder) -> Result<Vec<usize>> {
    match mode {
        ExecutionOrder::Insertion => Ok((0..flowsheet.instances.len()).collect()),
        ExecutionOrder::Topological => topological(flowsheet),
    }
}

/// The class's post-order DFS, over the instances in declaration order.
///
/// **The edges are the connections between instances, less the declared recycles.** A connection
/// with a bare name on either side is a feed or a product — a boundary — and contributes no
/// instance-to-instance edge, which is the same reading `check` takes.
fn topological(flowsheet: &Flowsheet) -> Result<Vec<usize>> {
    let n = flowsheet.instances.len();
    let index: HashMap<&str, usize> = flowsheet
        .instances
        .iter()
        .enumerate()
        .map(|(i, instance)| (instance.id.as_str(), i))
        .collect();

    // The recycles, as (from, to) instance pairs, so the edge they name can be skipped.
    let mut torn: Vec<(String, String)> = Vec::new();
    for recycle in &flowsheet.recycles {
        if let (Some((from, _)), Some((to, _))) =
            (instance_of(&recycle.from), instance_of(&recycle.to))
        {
            torn.push((from.to_string(), to.to_string()));
        }
    }

    // The outgoing edges, in the order the connections are declared - which is what makes this
    // walk reproducible. The class iterates `node.getOutgoingEdges()`, so its tie-breaking is
    // the order the edges were added too.
    let mut outgoing: Vec<Vec<usize>> = vec![Vec::new(); n];
    for connection in &flowsheet.connections {
        let (Some((from, _)), Some((to, _))) =
            (instance_of(&connection.from), instance_of(&connection.to))
        else {
            continue;
        };
        if torn.iter().any(|(a, b)| a == from && b == to) {
            continue;
        }
        let (Some(&source), Some(&target)) = (index.get(from), index.get(to)) else {
            return Err(AzothError::invalid_input(
                "connections",
                format!(
                    "`{}` names an instance the flowsheet does not declare",
                    if index.contains_key(from) { to } else { from }
                ),
            ));
        };
        outgoing[source].push(target);
    }

    let mut state = vec![Visit::Unvisited; n];
    let mut stack: Vec<usize> = Vec::with_capacity(n);
    for node in 0..n {
        if state[node] == Visit::Unvisited {
            visit(node, &outgoing, &mut state, &mut stack)?;
        }
    }
    // The post-order stack is reversed for topological order, which is what the class does.
    stack.reverse();
    Ok(stack)
}

/// What the walk has done to a node.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Visit {
    /// Not reached.
    Unvisited,
    /// On the walk's own path, so an edge back to it closes a cycle.
    OnPath,
    /// Finished, and already on the stack.
    Done,
}

/// `topologicalSortDFS`, with the cycle it cannot tear refused rather than looped on.
///
/// **The class's version cannot terminate wrongly because `analyzeCycles` has already marked the
/// back edges.** azoth's recycles are *declared* rather than detected, so an undeclared cycle is
/// a flowsheet the checker already refuses (`UnrecycledLoop`) - and this walk refuses it too
/// rather than recursing forever, which is the same statement made where it would otherwise hang.
fn visit(
    node: usize,
    outgoing: &[Vec<usize>],
    state: &mut [Visit],
    stack: &mut Vec<usize>,
) -> Result<()> {
    state[node] = Visit::OnPath;
    for &next in &outgoing[node] {
        match state[next] {
            Visit::OnPath => {
                return Err(AzothError::invalid_input(
                    "connections",
                    "an instance cycle is not declared as a recycle, so the order it runs in is \
                     not defined",
                ));
            }
            Visit::Unvisited => visit(next, outgoing, state, stack)?,
            Visit::Done => {}
        }
    }
    state[node] = Visit::Done;
    stack.push(node);
    Ok(())
}

/// The instance an endpoint names, or `None` for a boundary.
///
/// `from` is a feed name or `instance.port`; `to` is a product name or `instance.port`. A bare
/// name is a boundary on both sides, and the dot is what tells them apart - `check`'s own
/// reading.
fn instance_of(endpoint: &str) -> Option<(&str, &str)> {
    endpoint.split_once('.')
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::flowsheet::{Connection, Instance, Recycle, named_input};

    fn instance(id: &str) -> Instance {
        Instance {
            id: id.to_string(),
            unit: "unit_ops.heater".to_string(),
            parameters: Default::default(),
        }
    }

    fn connection(from: &str, to: &str) -> Connection {
        Connection {
            from: from.to_string(),
            to: to.to_string(),
        }
    }

    /// Three units in a chain, declared in the order they run.
    fn chain() -> Flowsheet {
        Flowsheet {
            id: "chain".to_string(),
            name: "A chain".to_string(),
            inputs: vec![named_input("in")],
            products: vec!["out".to_string()],
            instances: vec![instance("a"), instance("b"), instance("c")],
            connections: vec![
                connection("in", "a.inlet"),
                connection("a.outlet", "b.inlet"),
                connection("b.outlet", "c.inlet"),
                connection("c.outlet", "out"),
            ],
            recycles: Vec::new(),
        }
    }

    /// **The default is insertion order**, which is the measurement `useGraphBasedExecution`'s
    /// `false` gives.
    #[test]
    fn insertion_order_is_the_declaration_order() {
        let flowsheet = chain();
        assert_eq!(
            execution_order(&flowsheet, ExecutionOrder::Insertion).expect("it orders"),
            vec![0, 1, 2]
        );
        assert_eq!(ExecutionOrder::default(), ExecutionOrder::Insertion);
    }

    /// **Declared backwards, the two orders disagree** - which is what makes the mode a real
    /// choice rather than two spellings of one answer.
    #[test]
    fn topological_order_is_not_the_declaration_order() {
        let mut flowsheet = chain();
        flowsheet.instances = vec![instance("c"), instance("b"), instance("a")];
        let order = execution_order(&flowsheet, ExecutionOrder::Topological).expect("it orders");
        let ids: Vec<&str> = order
            .iter()
            .map(|&i| flowsheet.instances[i].id.as_str())
            .collect();
        assert_eq!(
            ids,
            vec!["a", "b", "c"],
            "the producer runs before its consumer"
        );
        assert_eq!(
            execution_order(&flowsheet, ExecutionOrder::Insertion).expect("it orders"),
            vec![0, 1, 2],
            "and insertion order is still the declaration's"
        );
    }

    /// **A declared recycle is not an ordering edge.** The class marks it a back edge and skips
    /// it in the walk; a port that kept it would put `b` before `a` and run a different machine.
    #[test]
    fn a_declared_recycle_does_not_constrain_the_order() {
        let flowsheet = Flowsheet {
            id: "loop".to_string(),
            name: "A loop".to_string(),
            inputs: vec![named_input("in")],
            products: vec!["out".to_string()],
            instances: vec![instance("a"), instance("b")],
            connections: vec![
                connection("in", "a.inlet"),
                connection("a.outlet", "b.inlet"),
                connection("b.outlet", "a.inlet"),
                connection("b.outlet", "out"),
            ],
            recycles: vec![Recycle::new("recycle_1", "b.outlet", "a.inlet")],
        };
        assert_eq!(
            execution_order(&flowsheet, ExecutionOrder::Topological).expect("it orders"),
            vec![0, 1]
        );
        // Without the declaration the same connections are a cycle, and the walk refuses it
        // rather than looping - the same statement `check`'s `UnrecycledLoop` makes.
        let mut untorn = flowsheet.clone();
        untorn.recycles.clear();
        assert!(execution_order(&untorn, ExecutionOrder::Topological).is_err());
    }
}
