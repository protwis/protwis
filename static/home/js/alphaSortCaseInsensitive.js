function localeCompareSupportsLocales1() {
  var a = "foo";
  var b = "bar";
  try {
    a.localeCompare(b, "en-us", {sensitivity:"base"});
  } catch (e) {
    return false;
  }
  return true;
}

function localeCompareSupportsLocales2() {
  var a = "foo";
  var b = "bar";
  try {
    a.localeCompare(b, "en", {sensitivity:"base"});
  } catch (e) {
    return false;
  }
  return true;
}

function localeCompareSupportsLocales3() {
  var a = "foo";
  var b = "bar";
  try {
    a.localeCompare(b, "en");
  } catch (e) {
    return false;
  }
  return true;
}


if(localeCompareSupportsLocales1()) {
  function alphaSortCaseInsensitive(a,b) {
    return a.localeCompare(b, "en-us", {sensitivity:"base"});
  }
} else if (localeCompareSupportsLocales2()) {
  function alphaSortCaseInsensitive(a,b) {
    return a.localeCompare(b, "en", {sensitivity:"base"});
  }
} else if (localeCompareSupportsLocales3()) {
  function alphaSortCaseInsensitive(a,b) {
    return a.localeCompare(b, "en");
  }
} else {
  function alphaSortCaseInsensitive(a,b) {
    return a.localeCompare(b);
  }
}