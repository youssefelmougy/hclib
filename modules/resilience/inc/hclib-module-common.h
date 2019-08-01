/*
 * Copy of modules/common/inc/hclib-module-common.h
 */

#ifndef HCLIB_MODULE_COMMON_H
#define HCLIB_MODULE_COMMON_H

namespace hclib {

template<class pending_op>
void poll_on_pending(pending_op **addr_of_head,
        bool (*test_completion_callback)(void *),
        hclib::locale_t *locale_to_yield_to) {
    do {
        int pending_list_non_empty = 1;

        pending_op *prev = NULL;
        pending_op *op = *addr_of_head;

        assert(op != NULL);

        while (op) {
            pending_op *next = op->next;

            const bool complete = test_completion_callback(op);

            if (complete) {
                // Remove from singly linked list
                if (prev == NULL) {
                    /*
                     * If previous is NULL, we *may* be looking at the front of
                     * the list. It is also possible that another thread in the
                     * meantime came along and added an entry to the front of
                     * this singly-linked wait list, in which case we need to
                     * ensure we update its next rather than updating the list
                     * head. We do this by first trying to automatically update
                     * the list head to be the next of wait_set, and if we fail
                     * then we know we have a new head whose next points to
                     * wait_set and which should be updated.
                     */
                    pending_op *old_head = __sync_val_compare_and_swap(
                            addr_of_head, op, op->next);
                    if (old_head != op) {
                        // Failed, someone else added a different head
                        assert(old_head->next == op);
                        old_head->next = op->next;
                        prev = old_head;
                    } else {
                        /*
                         * Success, new head is now wait_set->next. We want this
                         * polling task to exit if we just set the head to NULL.
                         * It is the responsibility of future async_when calls
                         * to restart it upon discovering a null head.
                         */
                        pending_list_non_empty = (op->next != NULL);
                    }
                } else {
                    /*
                     * If previous is non-null, we just adjust its next link to
                     * jump over the current node.
                     */
                    assert(prev->next == op);
                    prev->next = op->next;
                }

                if (op->serialized) {
                    op->data->deserialize(op->serialized);
                    //We do not want to delete data since the deserialize routine might
                    //directy use data pointer instead of copying it for efficiency
                    op->serialized->data = nullptr;
                    delete op->serialized;
                }
                else {
                    //Isend do not need to deserialize data once the operation is completed
                    //assert(false);
                }
                if (op->prom) {
                    hclib_promise_put(op->prom, op->data);
                } else {
                    //Allow the completion notification promise to be not set
                    //assert(false);
                }
                free(op);
            } else {
                prev = op;
            }

            op = next;
        }

        if (pending_list_non_empty) {
            hclib::yield_at(locale_to_yield_to);
        } else {
            // Empty list
            break;
        }
    } while (true);
}

template<class pending_op>
void append_to_pending(pending_op *op, pending_op **addr_of_head,
        bool (*test_completion_callback)(void *),
        hclib::locale_t *locale_to_yield_to) {
    pending_op *pending = *addr_of_head;
    op->next = pending;

    pending_op *old_head;
    while (1) {
        old_head = __sync_val_compare_and_swap(addr_of_head, op->next, op);
        if (old_head != op->next) {
            op->next = old_head;
        } else {
            break;
        }
    }

    if (old_head == NULL) {
        hclib::async_at([addr_of_head, test_completion_callback, locale_to_yield_to] {
            hclib::poll_on_pending<pending_op>(addr_of_head,
                test_completion_callback, locale_to_yield_to);
        }, locale_to_yield_to);
    }
}

}

#endif
