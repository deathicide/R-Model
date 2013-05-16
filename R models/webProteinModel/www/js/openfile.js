$(document).on("click", "button.open", function(evt) {

})

var openBinding = new Shiny.InputBinding();
$.extend(openBinding,{
  find: function(scope){
    return $(scope).find(".openFile");
  },
  subscribe: function(el, callback) {
    $(el).on("change.openBinding", function(e) {
      callback();
    });
  },
  unsubscribe: function(el) {
    $(el).off(".openBinding");
  }
});

Shiny.inputBindings.register(openBinding);