#include "graph.h"

/* create a new graph */
graph* new_graph()
{
     graph* g = ckallocz(sizeof(graph));
    g->node_hash = new_hashtable(8);
    return g;
}

/* create a new node for the graph. */
node* new_node(const char* const name, void* const val)
{
    node* n = ckallocz(sizeof(node));
    n->name = ckallocz(strlen(name) + 1);
    memcpy(n->name, name, strlen(name));
    n->val  = val;

    n->component = -1;    
    n->visited   = FALSE;    

    /* no edges initially */
    n->num_edges = 0;
    n->edges = NULL;

    return n;    
}

/* return TRUE if an edge exists from n1 to n2 */
static int edge_exists(const node* const n1, const node* const n2)
{
    pre(n1 != NULL);
    pre(n2 != NULL);
    pre(n1 != n2);

    int i;
    for(i = 0; i < n1->num_edges; i++){
        forceassert(n1->edges[i]->n1 == n1 || n1->edges[i]->n2 == n1);
        if(n1->edges[i]->n1 == n2 || n1->edges[i]->n2 == n2){
            return TRUE;
        }
    }

    return FALSE;
}

/* make an edge from n1 to n2 */
static void make_edge(graph* const g, 
                      node* const n1, 
                      node* const n2, 
                      const int weight)
{
    /* no self edges */
    pre(n1 != n2);
    pre(n1 != NULL);
    pre(n2 != NULL);

    /* does an edge exist from n1 to n2 */
    if(!edge_exists(n1, n2)){
        /* create the edge. */
        edge* e    = ckalloc(sizeof(edge));
        e->next    = NULL;

        e->n1   = n1;
        e->n2   = n2;
        e->weight = weight;
        
        /* add it both of the nodes edge list */
        n1->num_edges++;
        n1->edges = ckrealloc(n1->edges, n1->num_edges * sizeof(edge*));
        n1->edges[n1->num_edges - 1] = e;

        n2->num_edges++;
        n2->edges = ckrealloc(n2->edges, n2->num_edges * sizeof(edge*));
        n2->edges[n2->num_edges - 1] = e;

        sladdhead(&g->edge_list, e);
    }
}

/* add this node to the graph */
void add_node(graph* g, node* const n)
{
    pre(g != NULL);
    pre(n != NULL);

    // add this node to the hashtable of nodes.
    bin* bin = add_hashtable(g->node_hash, n->name, strlen(n->name), n);
    ckfree(bin->name);
    *(&bin->name) = n->name;
 
    // form the edges connecting this node to the other nodes in the graph. 
    evidence* evdnc1 = n->val;
    evidence* evdnc2 = NULL;
    node* iter = g->node_list;
    while(iter){
        evdnc2 = iter->val;
        forceassert(evdnc2->b1 <= evdnc1->b1);

        if((evdnc1->type == PAIRED_READ) && (evdnc2->type == PAIRED_READ)){
            // do these breakpoints overlap and support the same kind of
            // variant?
            if((evdnc2->b1 < evdnc1->b2) && 
               (evdnc1->variantclass == evdnc2->variantclass)){
                // if these support the same variant, then the composite
                // breakpoint would be [b1,b2]. Using those breakpoints, do both
                // of them satisfy the restriction for insert lengths for their
                // read groups?
                int32_t b1 = MAX(evdnc1->b1, evdnc2->b1);
                int32_t b2 = MIN(evdnc1->b2, evdnc2->b2);
                
                assert(evdnc1->aln1->start < b1);
                assert(evdnc2->aln1->start < b1);
                int32_t d1 = b1 - evdnc1->aln1->start 
                            + ((readseg*)sllast(evdnc1->aln3))->end-b2;
                int32_t d2 = b1 - evdnc2->aln1->start
                            + ((readseg*)sllast(evdnc2->aln3))->end-b2;
                if((d1 < evdnc1->max) && (d2 < evdnc2->max)){
                    make_edge(g, iter, n, 1);
                }
            }
        }else if((evdnc1->type == SPLIT_READ) && (evdnc2->type == SPLIT_READ)){
            if((evdnc1->variantclass == evdnc2->variantclass) && 
               (evdnc1->b1 == evdnc2->b1) && 
               (evdnc1->b2 == evdnc2->b2)){
                make_edge(g, iter, n, 1);
            }
        }else if((evdnc1->type == PAIRED_READ) && (evdnc2->type == SPLIT_READ)){
//            int32_t d1 = evdnc2->b1 - evdnc1->aln1->start
//                       + ((readseg*)sllast(evdnc1->aln3))->end - evdnc2->b2;
//            if(d1 < evdnc1->max){
//                make_edge(g, iter, n, 1);
//            }
        }else if((evdnc1->type == SPLIT_READ) && (evdnc2->type == PAIRED_READ)){ 
//            int32_t d2 = evdnc1->b1 - evdnc2->aln1->start
//                       + ((readseg*)sllast(evdnc2->aln3))->end - evdnc1->b2;
//            if(d2 < evdnc2->max){
//                make_edge(g, iter, n, 1);
//            }
        }else{
            fatalf("Unknown overlap of evidence types: %d and %d", evdnc1->type, evdnc2->type);
        }

        iter = iter->next;    
    }

    /* now add this to the node list */ 
    sladdhead(&g->node_list, n);
}

static int sort_by_component(const void* const el1,
                             const void* const el2)
{
    node* a = *((node**)el1);
    node* b = *((node**)el2);

    return a->component - b->component;
}

/* sort the nodes */
void sort_nodes(graph* g)
{
    slsort(&g->node_list, sort_by_component);
}

static void print_bin(bin* b)
{
    pre(b != NULL);

    printf("\"%s\";\n", b->name);
    node* n = b->val;

    int i;
    for(i = 0; i < n->num_edges; i++){
        printf("\"%s\" -- \"%s\"\n", 
              n->edges[i]->n1->name, n->edges[i]->n2->name);
    }
}

/* print a graph in the GRAPHVIZ format */
void print_graph(const graph* const g)
{
    pre(g != NULL);

    printf("graph g {\n");
    printf("node [shape = circle];\n");
    func_hashtable(g->node_hash, print_bin);
    printf("}\n");
}

static void depthfirst(node* const root, const int identifier)
{
    int i;
    node* next;
    root->visited = TRUE;
    assert(root->component == -1);
    root->component = identifier;

    for(i = 0; i < root->num_edges; i++){
        next = root->edges[i]->n1 == root 
             ? root->edges[i]->n2 : root->edges[i]->n1;

        if(next->visited == FALSE){
            depthfirst(next, identifier);   
        }
    }
}

/* find and mark the connected comonents in this graph. Return the number of
 * connected components in the graph */
int find_connected_components(graph* const graph)
{
    int componentId = 0;

    node* iter;
    for(iter = graph->node_list; iter; iter = iter->next){
        if(iter->visited == FALSE){
            /*mark all the nodes that you can reach from this node*/
            depthfirst(iter, componentId);
            componentId++;
        }
    }

   return componentId;
}

/* free the resources held by this graph */
void free_graph(graph** pg)
{
    graph* g = *pg;
    
    free_hashtable(&g->node_hash);

    node* iter = g->node_list;
    for(; iter; iter = iter->next){
        ckfree(iter->edges);
    }

    slfreelist(&g->node_list);  
    slfreelist(&g->edge_list);

    ckfree(g);
}

