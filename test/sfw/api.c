/*
 * api.c - API wrappers
 */

#include "ptl_test.h"

int test_ptl_init(struct node_info *info)
{
	libtest_init();
	info->err = PtlInit();

	return info->ret != info->err;
}

void test_ptl_fini(struct node_info *info)
{
	PtlFini();
	libtest_fini();
}

int test_ptl_ni_init(struct node_info *info)
{
	int ret;

	ret = PtlNIInit(info->iface, info->ni_opt, info->pid,
				      info->desired_ptr, info->actual_ptr,
				      info->ptr);

	info->err = ret;

	if (ret == PTL_OK && info->ptr == &info->ni_handle) {
		if (info->next_ni >= STACK_SIZE) {
			printf("NI stack overflow\n");
			return PTL_FAIL;
		}
		info->ni_stack[info->next_ni++] = info->ni_handle;
	}

	if (info->ni_opt & PTL_NI_LOGICAL) {
		PtlSetMap(info->ni_handle, info->map_size, info->desired_map_ptr);
	}

	return info->ret != ret;
}

int test_ptl_ni_fini(struct node_info *info)
{
	int ret;

	ret = PtlNIFini(info->ni_handle);

	info->err = ret;

	if (ret == PTL_OK) {
		if (info->next_ni == 0) {
			printf("NI stack underflow\n");
			return PTL_FAIL;
		}
		info->next_ni--;
	}

	return info->ret != ret;
}

int test_ptl_ni_status(struct node_info *info)
{
	int ret;

	ret = PtlNIStatus(info->ni_handle, info->reg, info->ptr);

	info->err = ret;

	return info->ret != ret;
}

int test_ptl_ni_handle(struct node_info *info)
{
	int ret;

	ret = PtlNIHandle(info->handle, info->ptr);

	info->err = ret;

	return info->ret != ret;
}

int test_ptl_handle_is_eq(struct node_info *info)
{
	int ret;

	ret = PtlHandleIsEqual(info->handle1, info->handle2);

	info->err = ret ? PTL_OK : PTL_FAIL;

	return info->ret != info->err;
}

int test_ptl_get_uid(struct node_info *info)
{
	int ret;

	ret = PtlGetUid(info->ni_handle, info->ptr);

	info->err = ret;

	return info->ret != ret;
}

int test_ptl_get_id(struct node_info *info)
{
	int ret;

	ret = PtlGetId(info->ni_handle, info->ptr);

	info->err = ret;

	return info->ret != ret;
}

int test_ptl_pt_alloc(struct node_info *info)
{
	int ret;

	ret = PtlPTAlloc(info->ni_handle, info->pt_opt,
		         info->eq_handle, info->pt_index, info->ptr);

	info->err = ret;

	if (ret == PTL_OK && info->ptr == &info->pt_index) {
		if (info->next_pt >= STACK_SIZE) {
			printf("PT stack overflow\n");
			return PTL_FAIL;
		}
		info->pt_stack[info->next_pt++] = info->pt_index;
	}

	return info->ret != ret;
}

int test_ptl_pt_free(struct node_info *info)
{
	int ret;

	while ((ret = PtlPTFree(info->ni_handle, info->pt_index)) ==
		PTL_PT_IN_USE)
		sched_yield();

	info->err = ret;

	if (ret == PTL_OK) {
		if (info->next_pt == 0) {
			printf("PT stack underflow\n");
			return PTL_FAIL;
		}
		info->next_pt--;
	}

	return info->ret != ret;
}

int test_ptl_pt_disable(struct node_info *info)
{
	int ret;

	ret = PtlPTDisable(info->ni_handle, info->pt_index);

	info->err = ret;

	return info->ret != ret;
}

int test_ptl_pt_enable(struct node_info *info)
{
	int ret;

	ret = PtlPTEnable(info->ni_handle, info->pt_index);

	info->err = ret;

	return info->ret != ret;
}

int test_ptl_eq_alloc(struct node_info *info)
{
	int ret;

	ret = PtlEQAlloc(info->ni_handle, info->eq_count, info->ptr);

	info->err = ret;

	if (ret == PTL_OK && info->ptr == &info->eq_handle) {
		if (info->next_eq >= STACK_SIZE) {
			printf("EQ stack overflow\n");
			return PTL_FAIL;
		}
		info->eq_stack[info->next_eq++] = info->eq_handle;
	}

	return info->ret != ret;
}

int test_ptl_eq_free(struct node_info *info)
{
	int ret;

	ret = PtlEQFree(info->eq_handle);

	info->err = ret;

	if (ret == PTL_OK) {
		if (info->next_eq == 0) {
			printf("EQ stack underflow\n");
			return PTL_FAIL;
		}
		info->next_eq--;
	}

	return info->ret != ret;
}

int test_ptl_eq_get(struct node_info *info)
{
	int ret;

	ret = PtlEQGet(info->eq_handle, info->ptr);

	info->err = ret;

	return info->ret != ret;
}

int test_ptl_eq_wait(struct node_info *info)
{
	int ret;

	ret = PtlEQWait(info->eq_handle, info->ptr);

	info->err = ret;

	return info->ret != ret;
}

int test_ptl_eq_poll(struct node_info *info)
{
	int ret;

	ret = PtlEQPoll(&info->eq_handle, info->eq_size,
			info->timeout, info->ptr, info->which_ptr);

	info->err = ret;

	return info->ret != ret;
}

int test_ptl_ct_alloc(struct node_info *info)
{
	int ret;

	ret = PtlCTAlloc(info->ni_handle, info->ptr);

	info->err = ret;

	if (ret == PTL_OK && info->ptr == &info->ct_handle) {
		if (info->next_ct >= STACK_SIZE) {
			printf("CT stack overflow\n");
			return PTL_FAIL;
		}
		info->ct_stack[info->next_ct++] = info->ct_handle;
	}

	return info->ret != ret;
}

int test_ptl_ct_free(struct node_info *info)
{
	int ret;

	ret = PtlCTFree(info->ct_handle);

	info->err = ret;

	if (ret == PTL_OK) {
		if (info->next_ct == 0) {
			printf("CT stack underflow\n");
			return PTL_FAIL;
		}
		info->next_ct--;
	}

	return info->ret != ret;
}

int test_ptl_ct_get(struct node_info *info)
{
	int ret;

	ret = PtlCTGet(info->ct_handle, info->ptr);

	info->err = ret;

	return info->ret != ret;
}

int test_ptl_ct_wait(struct node_info *info)
{
	int ret;

	ret = PtlCTWait(info->ct_handle, info->ct_test, info->ptr);

	info->err = ret;

	return info->ret != ret;
}

int test_ptl_ct_poll(struct node_info *info)
{
	int ret;

	ret = PtlCTPoll(&info->ct_handle, &info->ct_test,
			info->ct_size, info->timeout, info->ptr,
			info->which_ptr);

	info->err = ret;

	return info->ret != ret;
}

int test_ptl_ct_set(struct node_info *info)
{
	int ret;

	ret = PtlCTSet(info->ct_handle, info->ct_event);

	info->err = ret;

	return info->ret != ret;
}

int test_ptl_ct_inc(struct node_info *info)
{
	int ret;

	ret = PtlCTInc(info->ct_handle, info->ct_event);

	info->err = ret;

	return info->ret != ret;
}

int test_ptl_ct_cancel_trig(struct node_info *info)
{
	int ret;

	ret = PtlCTCancelTriggered(info->ct_handle);

	info->err = ret;

	return info->ret != ret;
}

int test_ptl_md_bind(struct node_info *info)
{
	int ret;

	info->md.ct_handle = info->ct_handle;
	info->md.eq_handle = info->eq_handle;

	if (debug) {
		printf("test_ptl_md_bind - start(%p), length(%d)\n",
			info->md.start, (int)info->md.length);
	}
	ret = PtlMDBind(info->ni_handle, &info->md, info->ptr);

	info->err = ret;

	if (ret == PTL_OK && info->ptr == &info->md_handle) {
		if (info->next_md >= STACK_SIZE) {
			printf("MD stack overflow\n");
			return PTL_FAIL;
		}
		info->md_stack[info->next_md++] = info->md_handle;
	}

	return info->ret != ret;
}

int test_ptl_md_release(struct node_info *info)
{
	int ret;
	int i;

	for (i = 0; i < 100; i++) {
		ret = PtlMDRelease(info->md_handle);
		if (ret == PTL_ARG_INVALID)
			usleep(10000);
		else
			break;
	}

	if (ret == PTL_OK) {
		if (info->next_md == 0) {
			printf("MD stack underflow\n");
			return PTL_FAIL;
		}
		info->next_md--;
	}

	return info->ret != ret;
}

int test_ptl_le_append(struct node_info *info)
{
	int ret;

	info->le.ct_handle = info->ct_handle;
	info->le.uid = info->uid;

	ret = PtlLEAppend(info->ni_handle, info->pt_index,
					&info->le, info->list, info->user_ptr,
					info->ptr);

	info->err = ret;

	if (ret == PTL_OK && info->ptr == &info->le_handle) {
		if (info->next_le >= STACK_SIZE) {
			printf("LE stack overflow\n");
			return PTL_FAIL;
		}
		info->le_stack[info->next_le++] = info->le_handle;
	}

	return info->ret != ret;
}

int test_ptl_le_unlink(struct node_info *info)
{
	int ret;

	while ((ret = PtlLEUnlink(info->le_handle)) == PTL_IN_USE)
		sched_yield();

	info->err = ret;

	if (ret == PTL_OK) {
		if (info->next_le == 0) {
			printf("LE stack underflow\n");
			return PTL_FAIL;
		}
		info->next_le--;
	}

	return info->ret != ret;
}

int test_ptl_le_search(struct node_info *info)
{
	int ret;

	info->le.ct_handle = info->ct_handle;

	ret = PtlLESearch(info->ni_handle, info->pt_index,
			  &info->le, info->search_op, info->user_ptr);

	info->err = ret;

	return info->ret != ret;
}

int test_ptl_me_append(struct node_info *info)
{
	int ret;

	info->me.ct_handle = info->ct_handle;
	info->me.uid = info->uid;

	ret = PtlMEAppend(info->ni_handle, info->pt_index,
					&info->me, info->list, info->user_ptr,
					info->ptr);

	info->err = ret;

	if (ret == PTL_OK && info->ptr == &info->me_handle) {
		if (info->next_me >= STACK_SIZE) {
			printf("ME stack overflow\n");
			return PTL_FAIL;
		}
		info->me_stack[info->next_me++] = info->me_handle;
	}

	return info->ret != ret;
}

int test_ptl_me_unlink(struct node_info *info)
{
	int ret;

	while ((ret = PtlMEUnlink(info->me_handle)) == PTL_IN_USE)
		sched_yield();

	info->err = ret;

	if (ret == PTL_OK) {
		if (info->next_me == 0) {
			printf("ME stack underflow\n");
			return PTL_FAIL;
		}
		info->next_me--;
	}

	return info->ret != ret;
}

int test_ptl_me_search(struct node_info *info)
{
	int ret;

	info->me.ct_handle = info->ct_handle;

	ret = PtlMESearch(info->ni_handle, info->pt_index,
			  &info->me, info->search_op, info->user_ptr);

	info->err = ret;

	return info->ret != ret;
}

int test_ptl_put(struct node_info *info)
{
	int ret;

	ret = PtlPut(info->md_handle, info->loc_offset,
		     info->length, info->ack_req, info->target_id,
		     info->pt_index, info->match,
		     info->rem_offset, info->user_ptr,
		     info->hdr_data);

	info->err = ret;

	return info->ret != ret;

}

int test_ptl_get(struct node_info *info)
{
	int ret;

	ret = PtlGet(info->md_handle, info->loc_offset,
		     info->length, info->target_id, info->pt_index,
		     info->match, info->rem_offset, info->user_ptr);

	info->err = ret;

	return info->ret != ret;
}

int test_ptl_atomic(struct node_info *info)
{
	int ret;

	ret = PtlAtomic(info->md_handle, info->loc_offset,
			info->length, info->ack_req, info->target_id,
			info->pt_index, info->match,
			info->rem_offset, info->user_ptr,
			info->hdr_data, info->atom_op, info->type);

	info->err = ret;

	return info->ret != ret;

}

int test_ptl_atomic_sync(struct node_info *info)
{
	int ret;

	ret = PtlAtomicSync();

	info->err = ret;

	return info->ret != ret;
}

int test_ptl_fetch_atomic(struct node_info *info)
{
	int ret;

	ret = PtlFetchAtomic(info->get_md_handle,
			     info->loc_get_offset,
			     info->put_md_handle,
			     info->loc_put_offset,
			     info->length, info->target_id,
			     info->pt_index, info->match,
			     info->rem_offset, info->user_ptr,
			     info->hdr_data, info->atom_op,
			     info->type);

	info->err = ret;

	return info->ret != ret;
}

int test_ptl_swap(struct node_info *info)
{
	int ret;

	ret = PtlSwap(info->get_md_handle, info->loc_get_offset,
		      info->put_md_handle, info->loc_put_offset,
		      info->length, info->target_id, info->pt_index,
		      info->match, info->rem_offset,
		      info->user_ptr, info->hdr_data, info->ptr,
		      info->atom_op, info->type);

	info->err = ret;

	return info->ret != ret;
}

int test_ptl_trig_put(struct node_info *info)
{
	int ret;

	ret = PtlTriggeredPut(info->md_handle, info->loc_offset,
			      info->length, info->ack_req,
			      info->target_id, info->pt_index,
			      info->match, info->rem_offset,
			      info->user_ptr, info->hdr_data,
			      info->ct_handle, info->threshold);

	info->err = ret;

	return info->ret != ret;
}

int test_ptl_trig_get(struct node_info *info)
{
	int ret;

	ret = PtlTriggeredGet(info->md_handle, info->loc_offset,
			      info->length, info->target_id,
			      info->pt_index, info->match,
			      info->rem_offset, info->user_ptr,
			      info->ct_handle, info->threshold);

	info->err = ret;

	return info->ret != ret;
}

int test_ptl_trig_atomic(struct node_info *info)
{
	int ret;

	ret = PtlTriggeredAtomic(info->md_handle, info->loc_offset,
				 info->length, info->ack_req,
				 info->target_id, info->pt_index, info->match,
				 info->rem_offset, info->user_ptr,
				 info->hdr_data, info->atom_op, info->type,
				 info->ct_handle, info->threshold);

	info->err = ret;

	return info->ret != ret;
}

int test_ptl_trig_fetch_atomic(struct node_info *info)
{
	int ret;

	ret = PtlTriggeredFetchAtomic(info->get_md_handle, info->loc_get_offset,
				      info->put_md_handle, info->loc_put_offset,
				      info->length, info->target_id,
				      info->pt_index, info->match,
				      info->rem_offset, info->user_ptr,
				      info->hdr_data, info->atom_op,
				      info->type, info->ct_handle,
				      info->threshold);

	info->err = ret;

	return info->ret != ret;
}

int test_ptl_trig_swap(struct node_info *info)
{
	int ret;

	ret = PtlTriggeredSwap(info->get_md_handle, info->loc_get_offset,
			       info->put_md_handle, info->loc_put_offset,
			       info->length, info->target_id, info->pt_index,
			       info->match, info->rem_offset,
			       info->user_ptr, info->hdr_data, info->ptr,
			       info->atom_op, info->type,
			       info->ct_handle, info->threshold);

	info->err = ret;

	return info->ret != ret;
}

int test_ptl_trig_ct_inc(struct node_info *info)
{
	int ret;

	ret = PtlTriggeredCTInc(info->ct_handle, info->ct_event,
				info->trig_ct_handle, info->threshold);

	info->err = ret;

	return info->ret != ret;
}

int test_ptl_trig_ct_set(struct node_info *info)
{
	int ret;

	ret = PtlTriggeredCTSet(info->ct_handle, info->ct_event,
				info->trig_ct_handle, info->threshold);

	info->err = ret;

	return info->ret != ret;
}

int test_ptl_start_bundle(struct node_info *info)
{
	int ret;

	ret = PtlStartBundle(info->ni_handle);

	info->err = ret;

	return info->ret != ret;
}

int test_ptl_end_bundle(struct node_info *info)
{
	int ret;

	ret = PtlEndBundle(info->ni_handle);

	info->err = ret;

	return info->ret != ret;
}

int test_ptl_set_map(struct node_info *info)
{
	int ret;

	ret = PtlSetMap(info->ni_handle, info->map_size, info->mapping);

	info->err = ret;

	return info->ret != ret;
}

int test_ptl_get_map(struct node_info *info)
{
	int ret;

	if (info->mapping)
		free(info->mapping);

	if (info->get_map_size)
		info->mapping = malloc(info->get_map_size * sizeof(ptl_process_t));
	else
		info->mapping = malloc(1);

	ret = PtlGetMap(info->ni_handle, info->get_map_size,
			info->mapping, &info->actual_map_size);

	info->err = ret;

	return info->ret != ret;
}
