#include "../headers/general.h"
#include "../headers/sequential.h"

#include <stdio.h>
#include <stdlib.h>

void track_group_init() {
    for (uint64_t i = 0; i < MAX_TRACKED_GROUPS; i++) {
        groups[i].count = 0;
    }
}

// -- Track a general object globally
void *track_generic(void *ptr, TrackType type) {
    if (!ptr || tracked_count >= MAX_TRACKED)
        return NULL;
    tracked[tracked_count++] = (TrackedObject){ptr, type};
    return ptr;
}

// -- Clear all globally tracked objects
void track_clear_all() {
    for (uint64_t i = 0; i < tracked_count; ++i) {
        switch (tracked[i].type) {
        case TRACK_MATRIX:
            m_free((matrix *)tracked[i].ptr);
            break;
        case TRACK_VECTOR:
            v_free((vector *)tracked[i].ptr);
            break;
        case TRACK_TENSOR:
            t_free((tensor *)tracked[i].ptr);
            break;
        }
        tracked[i].ptr = NULL;
    }
    tracked_count = 0;
}

// -- Track an object into a specific group
void *track_generic_group_in(void *ptr, uint64_t group_id, TrackType type) {
    if (!ptr || group_id >= MAX_TRACKED_GROUPS)
        return NULL;

    TrackedGroup *g = &groups[group_id];
    if (g->count < MAX_TRACKED_IN_GROUPS) {
        g->item[g->count++] = (TrackedObject){.ptr = ptr, .type = type};
        return ptr;
    }
    return NULL;
}

// -- Get a pointer to the object at index in group
void *track_generic_group_get(uint64_t group_id, uint64_t index) {
    if (group_id >= MAX_TRACKED_GROUPS)
        return NULL;

    TrackedGroup *g = &groups[group_id];
    if (index < g->count) {
        return g->item[index].ptr;
    }
    return NULL;
}

// -- Replace a tracked object in-place
void track_generic_group_replace(uint64_t group_id, uint64_t index, void *ptr,
                                 TrackType type) {
    if (!ptr || group_id >= MAX_TRACKED_GROUPS)
        return;

    TrackedGroup *g = &groups[group_id];
    if (index < g->count) {
        // Free previous object
        switch (g->item[index].type) {
        case TRACK_MATRIX:
            m_free((matrix *)g->item[index].ptr);
            break;
        case TRACK_VECTOR:
            v_free((vector *)g->item[index].ptr);
            break;
        case TRACK_TENSOR:
            t_free((tensor *)g->item[index].ptr);
            break;
        }
        g->item[index].ptr = ptr;
        g->item[index].type = type;
    }
}

// -- Clear and free all objects in a group
void track_group_clear(uint64_t group_id) {
    if (group_id >= MAX_TRACKED_GROUPS)
        return;

    TrackedGroup *g = &groups[group_id];
    for (uint64_t i = 0; i < g->count; i++) {
        switch (g->item[i].type) {
        case TRACK_MATRIX:
            m_free((matrix *)g->item[i].ptr);
            break;
        case TRACK_VECTOR:
            v_free((vector *)g->item[i].ptr);
            break;
        case TRACK_TENSOR:
            t_free((tensor *)g->item[i].ptr);
            break;
        }
        g->item[i].ptr = NULL;
        g->item[i].type = 0;
    }
    g->count = 0;
}

// -- Count number of tracked items in group
uint64_t track_group_count(uint64_t group_id) {
    if (group_id >= MAX_TRACKED_GROUPS)
        return 0;
    return groups[group_id].count;
}

uint64_t find_empty_group_id() {
    for (uint64_t i = 0; i < MAX_TRACKED_GROUPS; i++) {
        if (groups[i].count == 0)
            return i;
    }
    fprintf(stderr, "ERROR : no groups tracker disponible\n");
    return -1;
}
