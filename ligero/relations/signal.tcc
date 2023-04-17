#ifndef SIGNAL_TCC_
#define SIGNAL_TCC_

namespace ligero{
template<typename FieldT>
void signal<FieldT>::print(){
    if(this->is_primary_input)printf("a_%d ",this->index);
    else printf("b_%d ",this->index);
}
}
#endif
